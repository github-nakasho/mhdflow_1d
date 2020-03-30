#!/usr/bin/env/ python3

import numpy as np


class Ppm:
    def __init__(self, ix, order):
        order = 3
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))

    def mm(self, a, b):
        return 0.5 * (np.sign(a)+np.sign(b))*min(abs(a), abs(b))

    def make_rec(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        for m in range(8):
            for i in range(2, ix+2):
                dpp = V[m][i+2] - V[m][i+1]
                dp = V[m][i+1] - V[m][i]
                dm = V[m][i] - V[m][i-1]
                dmm = V[m][i-1] - V[m][i-2]
                d0bar = self.mm(0.5*(dp+dm), 2*self.mm(dm, dp))
                dpbar = self.mm(0.5*(dpp+dp), 2*self.mm(dp, dpp))
                dmbar = self.mm(0.5*(dm+dmm), 2*self.mm(dmm, dm))
                ddp = 0.5 * (V[m][i]+V[m][i+1]) - (dpbar-d0bar) / 6 - V[m][i]
                ddm = 0.5 * (V[m][i]+V[m][i-1]) + (dmbar-d0bar) / 6 - V[m][i]
                tmp1 = ddp
                tmp2 = ddm
                if tmp1 * tmp2 > 0:
                    ddp = 0
                    ddm = 0
                if abs(tmp1) >= 2 * abs(tmp2):
                    ddp = -2 * ddm
                elif abs(tmp2) >= 2 * abs(tmp1):
                    ddm = -2 * ddp                    
                Vl[m][i] = V[m][i] + ddp
                Vr[m][i-1] = V[m][i] + ddm
        return Vl, Vr
