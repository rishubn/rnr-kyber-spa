import time
import ctypes

from dataclasses import dataclass

import numpy as np
from picosdk.ps6000 import ps6000 as ps
from picosdk.functions import adc2mV, assert_pico_ok, mV2adc


ranges_mV = {
    10: 0,
    20: 1,
    50: 2,
    100: 3,
    200: 4,
    500: 5,
    1000: 6,
    2000: 7,
    5000: 8,
    10000: 9,
    20000: 10,
    50000: 11,
    100000: 12,
    200000: 13,
}

inv_ranges_mV = {v: k for k, v in ranges_mV.items()}


PS6000_RATIO_MODE = ps.PS6000_RATIO_MODE


@dataclass
class Picoscope6000:
    n_traces: int  # number of traces per acquisition
    # TODO trig delay
    range_signal: int  # mV, channel A
    post_trigger: int  # number of samples after the trigger
    pre_trigger: int = 0  # number of samples before trigger
    range_trigger: int = 5000  # mV, channel B
    trigger_threshold: int = 1600  # mV
    timebase: int = 3  # Set Picoscope's API guide.
    with_trigger: bool = False  # Save the trigger signal
    compensate_jitter: bool = (
        False  # Try to improve trace alignment with trigger-based picoscope timing information.
    )
    ratio_mode: int = PS6000_RATIO_MODE["PS6000_RATIO_MODE_NONE"]
    downsample_ratio: int = 0
    analogue_offset: float = 0.0

    def __enter__(self):
        self.setup_pico()
        return self

    def __exit__(self, type, value, tb):
        ps.ps6000CloseUnit(self.chandle)

    def setup_pico(self):
        # Create chandle and status ready for use
        self.chandle = ctypes.c_int16()
        self.status = {}
        ps.ps6000CloseUnit(self.chandle)

        ############## Step 1 #####################
        # Open the oscilloscope using ps6000OpenUnit.
        self.status["openunit"] = ps.ps6000OpenUnit(ctypes.byref(self.chandle), None)
        assert_pico_ok(self.status["openunit"])
        ch_a_range = ranges_mV[self.range_signal]
        self.maxvolt = ctypes.c_float(0.0)
        self.minvolt = ctypes.c_float(0.0)
        self.status["GetAnalogueOffset"] = ps.ps6000GetAnalogueOffset(
            self.chandle,
            ch_a_range,
            2,
            ctypes.byref(self.maxvolt),
            ctypes.byref(self.minvolt),
        )
        assert_pico_ok(self.status["GetAnalogueOffset"])
        assert self.analogue_offset < np.float32(
            self.maxvolt
        ) and self.analogue_offset > np.float32(self.minvolt)
        ############## Step 2 #####################
        # Select channel ranges and AC/DC coupling using ps6000SetChannel.

        self.status["setChA"] = ps.ps6000SetChannel(
            self.chandle, 0, 1, 0, ch_a_range, self.analogue_offset, 0
        )
        assert_pico_ok(self.status["setChA"])

        ch_b_range = ranges_mV[self.range_trigger]
        self.status["setChB"] = ps.ps6000SetChannel(
            self.chandle, 1, 1, 0, ch_b_range, 0, 0
        )
        assert_pico_ok(self.status["setChB"])

        self.status["SetExternalClock"] = ps.ps6000SetExternalClock(
            self.chandle, 2, 10000
        )
        assert_pico_ok(self.status["SetExternalClock"])

        ############## Step 3 #####################
        # Set the number of memory segments equal to or greater than the number of captures required using ps6000MemorySegments. Use ps6000SetNoOfCaptures before each run to specify the number of waveforms to capture.
        maxSamples = self.pre_trigger + self.post_trigger
        self.cmaxSamples = ctypes.c_int32(maxSamples)
        self.status["MemorySegments"] = ps.ps6000MemorySegments(
            self.chandle, self.n_traces, ctypes.byref(self.cmaxSamples)
        )
        assert_pico_ok(self.status["MemorySegments"])

        self.status["SetNoOfCaptures"] = ps.ps6000SetNoOfCaptures(
            self.chandle, self.n_traces
        )
        assert_pico_ok(self.status["SetNoOfCaptures"])

        ############## Step 4 #####################
        # Using ps6000GetTimebase, select timebases until the required nanoseconds per sample is located.
        timeIntervalns = ctypes.c_float()
        returnedMaxSamples = ctypes.c_int32()
        self.status["getTimebase2"] = ps.ps6000GetTimebase2(
            self.chandle,
            self.timebase,
            maxSamples,
            ctypes.byref(timeIntervalns),
            1,
            ctypes.byref(returnedMaxSamples),
            0,
        )
        assert_pico_ok(self.status["getTimebase2"])

        ############## Step 5 #####################
        # Use the trigger setup functions ps6000SetTriggerChannelConditions, ps6000SetTriggerChannelDirections and ps6000SetTriggerChannelProperties to set up the trigger if required.
        threshold_adc = int(
            mV2adc(self.trigger_threshold, ch_b_range, ctypes.c_int16(32512))
        )
        self.status["trigger"] = ps.ps6000SetSimpleTrigger(
            self.chandle, 1, 1, threshold_adc, 2, 0, 100000
        )
        assert_pico_ok(self.status["trigger"])

        ############## Step 8 #####################
        # Use ps6000SetDataBufferBulk to tell the driver where your memory buffers are. Call the function once for each channel/segment combination for which you require data. For greater efficiency with multiple captures, you could do this outside the loop after step 5.
        self.bufferAMax = []
        self.bufferAMin = []
        self.bufferBMax = []
        self.bufferBMin = []

        for i in range(self.n_traces):
            self.bufferAMax.append((ctypes.c_int16 * maxSamples)())
            self.bufferAMin.append((ctypes.c_int16 * maxSamples)())
            self.status["SetDataBuffersBulk"] = ps.ps6000SetDataBuffersBulk(
                self.chandle,
                0,
                ctypes.byref(self.bufferAMax[i]),
                ctypes.byref(self.bufferAMin[i]),
                maxSamples,
                i,
                self.ratio_mode,
            )
            assert_pico_ok(self.status["SetDataBuffersBulk"])

            if self.with_trigger:
                self.bufferBMax.append((ctypes.c_int16 * maxSamples)())
                self.bufferBMin.append((ctypes.c_int16 * maxSamples)())
                self.status["SetDataBuffersBulk"] = ps.ps6000SetDataBuffersBulk(
                    self.chandle,
                    1,
                    ctypes.byref(self.bufferBMax[i]),
                    ctypes.byref(self.bufferBMin[i]),
                    maxSamples,
                    i,
                    self.ratio_mode,
                )
            assert_pico_ok(self.status["SetDataBuffersBulk"])

    def arm_scope(self):
        ############## Step 6 #####################
        # Start the oscilloscope running using ps6000RunBlock.
        self.status["runBlock"] = ps.ps6000RunBlock(
            self.chandle,
            self.pre_trigger,
            self.post_trigger,
            self.timebase,
            0,
            None,
            0,
            None,
            None,
        )
        assert_pico_ok(self.status["runBlock"])

    def get_traces(self):
        ############## Step 7 #####################
        # Wait until the oscilloscope is ready using the ps6000BlockReady callback.
        ready = ctypes.c_int16(0)
        check = ctypes.c_int16(0)
        while ready.value == check.value:
            time.sleep(0.01)
            self.status["isReady"] = ps.ps6000IsReady(self.chandle, ctypes.byref(ready))

        ############## Step 9 #####################
        # Transfer the blocks of data from the oscilloscope using ps6000GetValuesBulk.
        overflow = (ctypes.c_int64 * self.n_traces)()
        self.status["GetValuesBulk"] = ps.ps6000GetValuesBulk(
            self.chandle,
            ctypes.byref(self.cmaxSamples),
            0,
            self.n_traces - 1,
            self.downsample_ratio,
            self.ratio_mode,
            ctypes.byref(overflow),
        )
        assert_pico_ok(self.status["GetValuesBulk"])

        ############# Step 10 ######################
        # Retrieve the time offset for each data segment using ps6000GetValuesTriggerTimeOffsetBulk64.
        times = (ctypes.c_int64 * self.n_traces)()
        timeUnits = (ctypes.c_int64 * self.n_traces)()
        self.status["GetValuesTriggerTimeOffsetBulk"] = (
            ps.ps6000GetValuesTriggerTimeOffsetBulk64(
                self.chandle,
                ctypes.byref(times),
                ctypes.byref(timeUnits),
                0,
                self.n_traces - 1,
            )
        )
        assert_pico_ok(self.status["GetValuesTriggerTimeOffsetBulk"])
        trig_offsets = np.array(times)

        data = np.array(self.bufferAMax, dtype=np.int16)

        if self.compensate_jitter:
            assert self.timebase == 3  # Otherwise trace alignment is destroyed.
            shifted = 0
            for i in range(0, self.n_traces):
                if (trig_offsets[i] // 1000) > 800:
                    data[i] = np.roll(data[i], 1)
                    shifted += 1

        if self.with_trigger:
            trigger = np.array(self.bufferBMax, dtype=np.int16)
            if self.compensate_jitter:
                for i in range(0, self.n_traces):
                    if (trig_offsets[i] // 1000) > 800:
                        trigger[i] = np.roll(trigger[i], 1)
            return data, trigger, np.array(overflow)
        else:
            return data, None, np.array(overflow)
