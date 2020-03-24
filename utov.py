#!/usr/bin/env/ python3

import numpy as np

def UtoV(U):
    V = np.zeros(U.shape)
    gm = 5.0 / 3.0
    gm1 = gm - 1
    oneorho = 1.0 / U[0]
    b2 = U[4] * U[4] + U[5] * U[5] + U[6] * U[6]
    r2v2 = U[1] * U[1] + U[2] * U[2] + U[3] * U[3]
    V[0] = U[0]
    V[1] = U[1] * oneorho
    V[2] = U[2] * oneorho
    V[3] = U[3] * oneorho
    V[4] = U[4]
    V[5] = U[5]
    V[6] = U[6]
    V[7] = gm1 * (U[7]-0.5*(r2v2*oneorho+b2))
    return V
