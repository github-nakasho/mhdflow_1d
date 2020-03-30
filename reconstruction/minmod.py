#!/usr/bin/env/ python3

import numpy as np


class Minmod:
    def __init__(self, ix):
        order = 2
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))

    def make_rec(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        for m in range(8):
            for i in range(1, ix+1):
                a = V[m][i] - V[m][i-1]
                b = V[m][i+1] - V[m][i]
                grad = np.sign(a) * max(0, min(abs(a), np.sign(a)*b))
                Vl[m][i] = V[m][i] + 0.5 * grad
                Vr[m][i-1] = V[m][i] - 0.5 * grad
        return Vl, Vr
