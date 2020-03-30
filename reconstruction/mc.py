#!/usr/bin/env/ python3

import numpy as np


class MC:
    def __init__(self, ix, order):
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
                c = 0.25 * (b + a)
                if a * b <= 0.0: 
                    grad = 0
                else:
                    if a >= 0:
                        grad = min(a, min(b, c))
                    else:
                        grad = max(a, max(b, c))
                Vl[m][i] = V[m][i] + grad
                Vr[m][i-1] = V[m][i] - grad
        return Vl, Vr
