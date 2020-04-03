#!/usr/bin/env/ python3

import numpy as np

from characteristic import Characteristic


class CENO:
    def __init__(self, ix, order):
        order = 3
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))
        self.a = np.array([1, 0.7, 1])
        self.char = Characteristic(ix, order)
        
    def mm(self, a, b):
        return 0.5 * (np.sign(a)+np.sign(b)) * min(abs(a), abs(b))

    def make_rec(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        for m in range(8):
            for i in range(2, ix+2):
                v = V[m][i]
                vp1 = V[m][i+1]
                vp2 = V[m][i+2]
                vm1 = V[m][i-1]
                vm2 = V[m][i-2]
                dvp = vp1 - v
                dvm = v - vm1
                dv = self.mm(dvp, dvm)
                vl = v + dv * 0.5
                vr = v + dv * (-0.5)
                d2p = 0.5 * (vp2-v)
                d20 = 0.5 * (vp1-vm1)
                d2m = 0.5 * (v-vm2)
                d3p = vp2 - 2 * vp1 + v
                d30 = vp1 - 2 * v + vm1
                d3m = v - 2 * vm1 + vm2
                ql = np.zeros(3)
                qr = np.zeros(3)
                dl = np.zeros(3)
                dr = np.zeros(3)
                ql[0] = vp1 + d2p * (-0.5) + 0.5 * d3p * (-0.5) ** 2
                ql[1] = v + d20 * 0.5 + 0.5 * d30 * 0.5 ** 2
                ql[2] = vm1 + d2m * 1.5 + 0.5 * d3m * (1.5) ** 2
                qr[0] = vp1 + d2p * (-1.5) + 0.5 * d3p * (-1.5) ** 2
                qr[1] = v + d20 * (-0.5) + 0.5 * d30 * (-0.5) ** 2
                qr[2] = vm1 + d2m * 0.5 + 0.5 * d3m * 0.5 ** 2
                dl = self.a * (ql-vl)
                dr = self.a * (qr-vr)
                Vl[m][i] = ql[np.argmin(abs(dl))]
                Vr[m][i-1] = qr[np.argmin(abs(dr))]
        return Vl, Vr
