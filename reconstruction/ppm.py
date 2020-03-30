#!/usr/bin/env/ python3

import numpy as np


class Ppm:
    def __init__(self, ix):
        order = 3
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))

    def phi(self, p, m):
        if p == 0 and m == 0:
            return 0
        else:
            return 2 * p * m / (p**2+m**2)

    def Ppm(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        k = 1 / 3
        for m in range(8):
            for i in range(2, ix+2):
                dup = V[m][i+1] - V[m][i]
                dum = V[m][i] - V[m][i-1]
                Vl[m][i] = V[m][i] + 0.25 * self.phi(dup, dum) * ((1-k)*dum+(1+k)*dup)
                Vr[m][i-1] = V[m][i] - 0.25 * self.phi(dup, dum) * ((1-k)*dup+(1+k)*dum)
        return Vl, Vr
