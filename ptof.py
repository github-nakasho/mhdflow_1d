#!/usr/bin/env/ python3

import numpy as np

def PtoF(ro, vx, vy, vz, bx, by, bz, pt, en):
    F = np.zeros((8, len(ro)))
    F[0] = ro * vx
    F[1] = ro * vx * vx + pt - bx * bx
    F[2] = ro * vy * vx - by * bx
    F[3] = ro * vz * vx - bz * bx
    F[4] = 0
    F[5] = by * vx - vy * bx
    F[6] = bz * vx - vz * bx
    F[7] = (en+pt) * vx - (vy*by+vz*bz) * bx
    return F
