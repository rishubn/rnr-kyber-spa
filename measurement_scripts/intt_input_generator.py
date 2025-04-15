from nttwrapper import INTT, INTTImpl
import numpy as np


class INTTInputGenerator:

    def __init__(self, intt: INTT, fixed=False):
        self.intt = intt
        self.fixed = fixed
        if fixed:
            self.fixed_msg = self.gen_msg()

    def gen_msg(self):
        match self.intt.impl:
            case INTTImpl.SIGNED | INTTImpl.WIDE_SIGNED:
                return self._gen_msg_s()
            case INTTImpl.UNSIGNED | INTTImpl.WIDE_UNSIGNED:
                return self._gen_msg_u()

    def _gen_msg_s(self):
        return np.random.randint(
            low=self.intt.domain.start,
            high=self.intt.domain.end,
            size=(32,),
            dtype=np.int16,
        )

    def _gen_msg_u(self):
        return np.random.randint(
            low=self.intt.domain.start,
            high=self.intt.domain.end,
            size=(32,),
            dtype=np.uint16,
        )
