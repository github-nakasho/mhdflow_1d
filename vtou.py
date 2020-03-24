#!/usr/bin/env/ python3

import numpy as np

def VtoU(V):
    U = np.zeros(V.shape)    
    gm = 5.0 / 3.0
    oneogm1 = 1.0 / (gm-1)
    v2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3]
    b2 = V[4] * V[4] + V[5] * V[5] + V[6] * V[6]
    U[0] = V[0]
    U[1] = V[0] * V[1]
    U[2] = V[0] * V[2]
    U[3] = V[0] * V[3]
    U[4] = V[4]
    U[5] = V[5]
    U[6] = V[6]
    U[7] = 0.5 * (V[0]*v2+b2) + oneogm1 * V[7]
    return U
