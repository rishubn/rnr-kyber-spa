import chipwhisperer as cw
from dataclasses import dataclass
import subprocess
import time
from chipwhisperer.capture import scopes


@dataclass
class CW308:
    scope: scopes.ScopeTypes = None
    baud: int = 52103
    output_len: int = 64

    def __enter__(self):
        self.setup_cw308()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.target.dis()
        if exc_type is not None:
            return False
        return True

    def setup_cw308(self):
        self.target = cw.target(self.scope, cw.targets.SimpleSerial)
        self.target.baud = self.baud
        self.target.output_len = self.output_len

    def program_target(self, binfile: str):
        cw.program_target(self.scope, cw.programmers.STM32FProgrammer, binfile)

    def write_bytes(self, op: str, msg: bytearray):
        assert len(msg) <= 64, "msg must be less than 64 bytes"
        self.target.simpleserial_write(op, msg)

    def read_bytes(self, op: str, output_len: int) -> bytearray:
        return self.target.simpleserial_read(op, output_len)

    def reset_target(self):
        self.scope.io.nrst = "low"
        time.sleep(0.25)
        self.scope.io.nrst = "high_z"
        time.sleep(0.25)
