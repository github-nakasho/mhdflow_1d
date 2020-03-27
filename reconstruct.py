#!/usr/bin/env/ python3

import numpy as np


class Reconstruct:
    def __init__(self, ixmax):
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))

    def Minmod(self, V, ix, order, dx, dt):
        Vl = self.Vl
        Vr = self.Vr
        for m in range(8):
            for i in range(1, ix+1):
                a = V[m][i] - V[m][i-1]
                b = V[m][i+1] - V[m][i]
                if a * b <= 0.0: 
                    grad = 0.0
                else:
                    if a > 0.0: 
                        grad = min(a, b)
                    else:
                        grad = max(a, b)
                Vl[m][i] = V[m][i] + 0.5 * grad
                Vr[m][i-1] = V[m][i] - 0.5 * grad
        return Vl, Vr

    def phi(self, p, m):
        if p == 0 and m == 0:
            return 0
        else:
            return 2 * p * m / (p**2+m**2)

    def Ppm(self, V, ix, order, dx, dt):
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