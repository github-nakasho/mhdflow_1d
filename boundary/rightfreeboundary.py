#!/usr/bin/env/ python3

class RightFreeBoundary:
    def __init__(self, ix, order):
        self.order = order
        self.ix = ix

    def set_u(self, U):
        order = self.order
        ixmax = self.ix + 2 * (order-1)
        for m in range(8):
            for i in range(self.ix+(order-1), ixmax):
                U[m][i] = U[m][i-1]
        return U

    def set_vr(self, Vl, Vr):
        order = self.order
        ix = self.ix
        for m in range(8):
            Vr[m][ix+(order-2)] = Vl[m][ix+(order-2)]
        return Vr
