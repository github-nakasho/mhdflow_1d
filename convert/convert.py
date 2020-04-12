#!/usr/bin/env/ python3

import numpy as np

from common import Const


class Convert:
    def __init__(self):
        self.gm = Const.GAMMA

    def PtoF(self, ro, vx, vy, vz, bx, by, bz, pt, en):
        F = np.zeros((8, len(ro)))
        F[0] = ro * vx
        F[1] = ro * vx * vx + pt - bx * bx
        F[2] = ro * vy * vx - by * bx
        F[3] = ro * vz * vx - bz * bx
        F[4] = (en+pt) * vx - (vy*by+vz*bz) * bx
        F[5] = 0
        F[6] = by * vx - vy * bx
        F[7] = bz * vx - vz * bx
        return F

    def UtoV(self, U):
        V = np.zeros(U.shape)
        gm = self.gm
        gm1 = gm - 1
        oneoro = 1 / U[0]
        b2 = U[5] * U[5] + U[6] * U[6] + U[7] * U[7]
        r2v2 = U[1] * U[1] + U[2] * U[2] + U[3] * U[3]
        V[0] = U[0]
        V[1] = U[1] * oneoro
        V[2] = U[2] * oneoro
        V[3] = U[3] * oneoro
        V[4] = gm1 * (U[4]-0.5*(r2v2*oneoro+b2))
        V[5] = U[5]
        V[6] = U[6]
        V[7] = U[7]
        return V

    def UVtoF(self, U, V, b2, vb):
        F = np.zeros(U.shape)
        F[0] = U[1]
        F[1] = U[1] * V[1] + (V[4]+0.5*b2) - V[5] * V[5]
        F[2] = U[2] * V[1] - V[6] * V[5]
        F[3] = U[3] * V[1] - V[7] * V[5]
        F[4] = (U[4]+V[4]+0.5*b2) * V[1] - vb * V[5]
        F[5] = 0.0
        F[6] = V[6] * V[1] - V[2] * V[5]
        F[7] = V[7] * V[1] - V[3] * V[5]
        return F

    def VtoU(self, V):
        U = np.zeros(V.shape)    
        gm = self.gm
        oneogm1 = 1 / (gm-1)
        v2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3]
        b2 = V[5] * V[5] + V[6] * V[6] + V[7] * V[7]
        U[0] = V[0]
        U[1] = V[0] * V[1]
        U[2] = V[0] * V[2]
        U[3] = V[0] * V[3]
        U[4] = 0.5 * (V[0]*v2+b2) + oneogm1 * V[4]
        U[5] = V[5]
        U[6] = V[6]
        U[7] = V[7]
        return U
