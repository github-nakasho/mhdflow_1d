#!/usr/bin/env/ python3

import numpy as np

def Minmod(V, ixmax, order, dx, dt):
    Vl = np.zeros((8, ixmax-(order-1)))
    Vr = np.zeros((8, ixmax-(order-1)))
    for m in range(8):
        for i in range(1, ixmax-(order-1)):
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
