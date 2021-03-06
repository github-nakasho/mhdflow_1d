#!/usr/bin/env/ python3

import numpy as np

from common import Const

def SetDt(V, minlength, cfl):
    gm = Const.GAMMA
    v2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3]
    b2 = V[5] * V[5] + V[6] * V[6] + V[7] * V[7]
    oneocfldt = max(np.sqrt(v2+(gm*V[4]+b2)/V[0])/minlength)
    dt = cfl/oneocfldt
    return dt
