#!/usr/bin/env/ python3

import numpy as np

def SetDt(V, minlength, cfl):
    gm = 5.0 / 3.0
    v2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3]
    b2 = V[4] * V[4] + V[5] * V[5] + V[6] * V[6]
    oneocfldt = max(np.sqrt(v2+(gm*V[7]+b2)/V[0])/minlength)
    dt = cfl/oneocfldt
    return dt
