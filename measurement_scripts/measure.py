from picoscope6000_sca import Picoscope6000
import chipwhisperer as cw
import matplotlib.pyplot as plt
from cw308 import CW308
from cwlite import CWLite
import numpy as np
from time import sleep
from argparse import ArgumentParser
import sys
from scalib.metrics import SNR
from tqdm import tqdm
import aes
from scalib.preprocessing import Quantizer
from traces import Trace
from intt_input_generator import INTTInputGenerator
import pickle
from pathlib import Path
from math import ceil
from loguru import logger
import json
import time
from enum import Enum
from nttwrapper import INTTImpl, INTT, KYBERQ, WIDE_S_KYBERQ, WIDE_U_KYBERQ
from numcodecs import Blosc
import zarr

bin_file = "/home/rnagpal/chipwhisperer/hardware/victims/firmware/simpleserial-ntt/simpleserial-ntt-CW308_STM32F3.hex"


def parse_args():
    argparse = ArgumentParser()
    argparse.add_argument("--program", action="store_true")
    argparse.add_argument("--ntraces", type=int)
    argparse.add_argument("--block", type=int)
    argparse.add_argument("--outpath")
    argparse.add_argument("--trace_length", type=int)
    argparse.add_argument("--fixed-msg", action="store_true")
    argparse.add_argument("--verify", action="store_true")
    argparse.add_argument("--debug", action="store_true")
    argparse.add_argument("--plot", action="store_true")
    argparse.add_argument(
        "--intt-type",
        type=lambda c: INTTImpl[c.upper()],
        choices=list(INTTImpl),
        required=True,
    )
    return argparse.parse_args()


# pre_trigger = 0
# post_trigger = 25000
pre_offset = 0
# ntraces = 1000
# block = 1000
# ns = pre_trigger + post_trigger
ntt_length = 230000
trace_length = 230000
max_adc_samples = 23000
# snr = SNR(nc=3329)


def measure_aes(ntraces, cw308, pbar):
    labels = np.zeros((ntraces, 1), dtype=np.uint16)
    cw308.output_len = 16
    for i in range(ntraces):
        msg, labels[i] = aes.gen_labels(rand=True)
        msg_byte = bytearray(msg.byteswap().tobytes())
        cw308.write_bytes("q", msg_byte)
        cw308.read_bytes("r").hex()
        pbar.update(1)
    return labels


def opcode(intt_type: INTTImpl):
    match intt_type:
        case INTTImpl.SIGNED:
            return "a"
        case INTTImpl.WIDE_SIGNED:
            return "b"
        case INTTImpl.UNSIGNED:
            return "c"
        case INTTImpl.WIDE_UNSIGNED:
            return "d"
        case _:
            raise NotImplementedError(f"didnt implement for {intt_type}")


def gen_labels(intt_out, intt_type):
    labels = np.array(intt_out)
    match intt_type:
        case INTTImpl.SIGNED:
            labels[labels < 0] += KYBERQ
        case INTTImpl.WIDE_SIGNED:
            labels[labels < 0] += WIDE_S_KYBERQ
        case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
            return labels
        case _:
            raise NotImplementedError(f"didnt implement for {intt_type}")
    return labels


def get_np_type(intt_type):
    match intt_type:
        case INTTImpl.SIGNED | INTTImpl.WIDE_SIGNED:
            return np.int16
        case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
            return np.uint16


def main():
    args = parse_args()
    ntraces = args.ntraces
    block = args.block
    pbar = tqdm(total=ntraces)
    compressor = Blosc(cname="lz4")
    outpath = Path(args.outpath)
    logger.add(f"{args.outpath}/trace.log", level="INFO")
    if args.debug:
        logger.add(f"{args.outpath}/debug.log", level="DEBUG")
    logger.info(f"Output dir {args.outpath}")
    with (
        CWLite() as scope,
        CW308(scope=scope.scope, baud=38400, output_len=16) as cw308,
        zarr.open(outpath / "trace_set.zarr", mode="w") as trace_set,
    ):

        if args.program:
            logger.info("Programming new binary...")
            cw308.program_target(bin_file)
        else:
            logger.info("Using existing binary...")
            cw308.reset_target()
        trace_length = args.trace_length
        if args.trace_length <= max_adc_samples:
            scope.samples = trace_length
            repeat_measure = 1
            logger.info(
                f"\nADC - {scope.samples} samples\ntrace length {trace_length}\n"
            )
        else:
            scope.samples = max_adc_samples
            repeat_measure = ceil(trace_length / max_adc_samples)
            logger.info(
                f"\nADC - {scope.samples} samples\ntrace length {trace_length}\nRepeating measurement {repeat_measure} times"
            )

        # traces = trace_set.zeros(
        #     "traces",
        #     shape=(ntraces, trace_length),
        #     chunks=(block, None),
        #     dtype=np.float64,
        #     compressor=compressor,
        # )
        labels = trace_set.zeros(
            "labels",
            shape=(ntraces, 256, 9),
            chunks=(block, 1, 1),
            dtype=np.uint16,
            compressor=`compressor`,
        )
        qtraces = trace_set.zeros(
            "qtraces",
            shape=(ntraces, trace_length),
            chunks=(block, None),
            dtype=np.int16,
            compressor=compressor,
        )
        pt = trace_set.zeros(
            "plaintext",
            shape=(ntraces, 32),
            chunks=(block, None),
            dtype=get_np_type(args.intt_type),
            compressor=compressor,
        )
        intt = INTT(args.intt_type)
        intt_input_generator = INTTInputGenerator(intt, args.fixed_msg)
        assert ntraces % block == 0
        logger.info("Chunk size {} Total chunks {}", block, ntraces // block)

        for chunk in range(ntraces // block):
            traces_chunk = np.zeros((block, trace_length))
            pt_chunk = np.zeros((block, 32), dtype=np.int16)
            labels_chunk = np.zeros((block, 256, 9), dtype=np.uint16)
            for n in range(block):
                if args.fixed_msg:
                    pt_chunk[n] = intt_input_generator.fixed_msg
                else:
                    pt_chunk[n] = intt_input_generator.gen_msg()
                intt_out = intt.compute_full(pt_chunk[n])
                labels_chunk[n] = gen_labels(intt_out, intt.impl)
                msg_byte = bytearray(pt_chunk[n].byteswap().tobytes())
                logger.debug(f"PT:\n{pt_chunk[n]}")
                logger.debug(f"INTT OUT:\n{intt_out}")
                logger.debug(f"LABELS:\n{labels_chunk}")
                resp = None
                for i in range(repeat_measure):
                    (
                        traces_chunk[n, i * scope.samples : (i + 1) * scope.samples],
                        resp,
                    ) = scope.capture_trace(
                        cw308,
                        (i * scope.samples) + 10,
                        msg_byte,
                        opcode(args.intt_type),
                    )
                if args.verify:
                    logger.debug(f"RESPONSE:\n{resp}")
                    assert (resp == intt_out[0:32, -1]).all()
                pbar.update(1)

            logger.debug(f"Saving chunk {chunk}")
            # assert (
            #     np.max(traces_chunk) < 0.5 and np.min(traces_chunk) > -0.5
            # ), f"max {np.max(traces_chunk)} min {np.min(traces_chunk)} {np.argmin(traces_chunk)}"
            # traces[block * chunk : block * (chunk + 1)] = traces_chunk
            pt[block * chunk : block * (chunk + 1)] = pt_chunk
            labels[block * chunk : block * (chunk + 1)] = labels_chunk
            qtraces[block * chunk : block * (chunk + 1)] = Quantizer.fit(
                traces_chunk
            ).quantize(traces_chunk)
        if args.plot:
            plt.figure()
            plt.plot(traces_chunk[0])
            plt.savefig(f"{args.outpath}/plot_{chunk}.png", dpi=300)
        logger.debug("Trace collection done, writing file...")
        with open(f"{args.outpath}/trace_stats.json", "w") as f:
            json.dump(
                {
                    "trace_length": trace_length,
                    "block": block,
                    "num_traces": ntraces,
                    "num_chunks": ntraces // block,
                    "trace_type": "fixed" if args.fixed_msg else "random",
                    "elapsed_time": time.strftime(
                        "%H:%M:%S", time.gmtime(pbar.format_dict["elapsed"])
                    ),
                },
                f,
                indent=4,
            )
        # Path(f"{args.outpath}/trace_done").touch()

        # for chunk in range(0, ntraces, block):

        # traces = scope.capture_traces(cw308, block, invntt.INVNTT(), pbar)
        # if quantizer is None:
        #     quantizer = Quantizer.fit(traces.waves)
        # traces.waves = quantizer.quantize(traces.waves)
        # pickle.dump(traces, f)
        # snr.fit_u(traces.waves, traces.labels)
    # plt.plot(traces.waves[0])
    # plt.savefig("out.png")
    # fig, ax = plt.subplots()
    # for i in range(4):
    #     s = snr.get_snr()[i]

    #     ax.plot(s)
    #     ax.set_yscale("log")
    # fig.savefig("snr.png")
    # with CW308() as cw308, Picoscope6000(
    #     n_traces=block,
    #     range_signal=500,
    #     post_trigger=post_trigger,
    #     pre_trigger=pre_trigger,
    #     range_trigger=500,
    #     trigger_threshold=160,
    #     timebase=5,
    #     with_trigger=True,
    # ) as picoscpe:

    #     for chunk in range(0, ntraces, block):
    #         # traces = np.zeros((block, ns), dtype=np.int16)
    #         if args.program:
    #             cw308.program_target(bin_file)
    #         else:
    #             cw308.reset_target()

    #         picoscpe.arm_scope()
    #         sleep(0.5)

    #         labels = measure_aes(block, cw308, pbar)
    #         # msg, l = verify.gen_labels(rand=True)
    #         # msg_byte = bytearray(msg.byteswap().tobytes())
    #         # cw308.write_bytes("p", msg_byte)
    #         # cw308.read_bytes("r")
    #         # labels[i] = l[2, 1] if l[2, 1] >= 0 else l[2, 1] + 3329
    #         # pbar.update(1)
    #         data, trigger, overflow = picoscpe.get_traces()
    #         # assert not overflow[pre_offset:].any()

    #         snr.fit_u(data, labels)
    # fig, ax = plt.subplots()

    # ax.set_ylim(-32512, 32512)
    # for i in range(1):
    #     ax.plot(
    #         np.arange(pre_trigger + post_trigger - pre_offset),
    #         data[i, pre_offset:],
    #         trigger[i, pre_offset:],
    #     )
    # np.savez_compressed("out.npz", traces=data)
    # fig.savefig("out.png")
    # print(data)
    # print(trigger)
    # print(overflow)


if __name__ == "__main__":
    main()
# with CW308() as cw308:
#    cw308.program_target(bin_file)
#    msg = bytearray()
#    for _ in range(32):
#        x = 2
#        msg.extend(x.to_bytes(2, 'big',signed=x < 0))
#    cw308.write_bytes('p', msg)
#    print(cw308.read_bytes('r'))
