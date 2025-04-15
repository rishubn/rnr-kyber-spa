from dataclasses import dataclass
import chipwhisperer as cw
import numpy as np
from traces import Trace
import time


@dataclass
class CWLite:
    clkfreq: int = 7.37e6

    def __enter__(self):
        self.setup_cwlite()
        return self

    def __exit__(self, exec_type, exec_value, tb):
        self.scope.dis()
        if exec_type is not None:
            return False
        return True

    @property
    def samples(self):
        return self.scope.adc.samples

    @samples.setter
    def samples(self, s):
        self.scope.adc.samples = s

    def setup_cwlite(self):
        self.scope = cw.scope()
        self.scope.default_setup()
        self.scope.adc.samples = 1000
        self.scope.clock.clkgen_freq = self.clkfreq

    # self.scope.gain.db = 28

    def arm_scope(self):
        self.scope.arm()

    def capture_trace(self, target, offset, msg_byte, opcode):
        self.scope.adc.offset = offset
        self.scope.arm()

        target.write_bytes(opcode, msg_byte)
        self.scope.capture()
        resp = np.frombuffer(target.read_bytes("r", 64), dtype=">i2")
        trace = self.scope.get_last_trace(as_int=False)
        return trace, resp

    def capture_traces(self, target, ntraces, model, pbar):
        waves = np.zeros((ntraces, self.scope.adc.samples))
        pt = np.zeros((ntraces, 32), dtype=np.int16)
        resp = np.zeros((ntraces, 32), dtype=np.int16)
        labels = np.zeros((ntraces, 4), dtype=np.uint16)
        for i in range(ntraces):
            self.scope.arm()
            msg = model.gen_msg(rand=True)
            msg_byte = bytearray(msg.byteswap().tobytes())
            target.write_bytes("p", msg_byte)
            self.scope.capture()
            pt[i] = msg
            resp = np.frombuffer(target.read_bytes("r", model.output_len), dtype=">i2")
            waves[i] = self.scope.get_last_trace(as_int=False)
            labels[i] = model.gen_labels(msg)
            pbar.update(1)
        return Trace(waves, pt, resp, labels, model.key)
