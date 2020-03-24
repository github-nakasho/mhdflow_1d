#!/usr/bin/env/ python3

import numpy as np

def Minmod(V, ixmax):
    Vl = np.zeros((8, ixmax))
    Vr = np.zeros((8, ixmax))
    for m in range(8):
        for i in range(1, ixmax-1):
            a = V[m][i] - V[m][i-1]
            b = V[m][i+1] - V[m][i]
            if a * b < 0.0: 
                grad = min(a, b)
            else:
                grad = max(a, b)
            Vl[m][i] = V[m][i] + 0.5 * grad
            Vr[m][i-1] = V[m][i] - 0.5 * grad
    return Vl, Vr
