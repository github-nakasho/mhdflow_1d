#!/usr/bin/env/ python3

class LeftFreeBoundary:
    def __init__(self, order):
        self.order = order

    def set_u(self, U):
        order = self.order
        for m in range(8):
            for i in range(order-1, 0, -1):
                U[m][i-1] = U[m][i]
        return U

    def set_vl(self, Vl, Vr):
        order = self.order
        for m in range(8):
            Vl[m][order-2] = Vr[m][order-2]
        return Vl
