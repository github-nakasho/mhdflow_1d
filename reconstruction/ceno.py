#!/usr/bin/env/ python3

import numpy as np


class CENO:
    def __init__(self, ix, order):
        order = 3
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))
        self.a = np.array([1, 0.7, 1])

    def mm(self, a, b):
        return 0.5 * (np.sign(a)+np.sign(b)) * min(abs(a), abs(b))

    def make_rec(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        for m in range(8):
            for i in range(2, ix+2):
                dvp = V[m][i+1]-V[m][i]
                dvm = V[m][i]-V[m][i-1]
                dv = self.mm(dvp, dvm)
                vl = V[m][i] + dv * 0.5
                vr = V[m][i] + dv * (-0.5)
                d2p = 0.5 * (V[m][i+2]-V[m][i])
                d20 = 0.5 * (V[m][i+1]-V[m][i-1])
                d2m = 0.5 * (V[m][i]-V[m][i-2])
                d3p = V[m][i+2] - 2 * V[m][i+1] + V[m][i]
                d30 = V[m][i+1] - 2 * V[m][i] + V[m][i-1]
                d3m = V[m][i] - 2 * V[m][i-1] + V[m][i-2]
                ql = np.zeros(3)
                qr = np.zeros(3)
                dl = np.zeros(3)
                dr = np.zeros(3)
                ql[0] = V[m][i+1] + d2p * (-0.5) + 0.5 * d3p * (-0.5) ** 2
                ql[1] = V[m][i] + d20 * 0.5 + 0.5 * d30 * 0.5 ** 2
                ql[2] = V[m][i-1] + d2m * 1.5 + 0.5 * d3m * (1.5) ** 2
                qr[0] = V[m][i+1] + d2p * (-1.5) + 0.5 * d3p * (-1.5) ** 2
                qr[1] = V[m][i] + d20 * (-0.5) + 0.5 * d30 * (-0.5) ** 2
                qr[2] = V[m][i-1] + d2m * 0.5 + 0.5 * d3m * 0.5 ** 2
                dl = self.a * (ql-vl)
                dr = self.a * (qr-vr)
                Vl[m][i] = ql[np.argmin(abs(dl))]
                Vr[m][i-1] = qr[np.argmin(abs(dr))]
        return Vl, Vr
