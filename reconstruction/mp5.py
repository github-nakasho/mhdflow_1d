#!/usr/bin/env/ python3

import numpy as np

from characteristic import Characteristic
from common import Const


class MP5:
    def __init__(self, ix, order):
        order = 5
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Wl = np.zeros((8, ixmax-1))
        self.Wr = np.zeros((8, ixmax-1))
        self.a = np.array([1, 0.7, 1])
        self.char = Characteristic(ix, order)

    def median(self, a, b, c):
        return a + self.mm(b-a, c-a)

    def mm(self, a, b):
        return 0.5 * (np.sign(a)+np.sign(b)) * min(abs(a), abs(b))
    
    def mm4(self, a, b, c, d):
        return 0.125 * (np.sign(a)+np.sign(b)) * abs((np.sign(a)+np.sign(b))*(np.sign(c)+np.sign(d))) * min(abs(a), abs(b), abs(c), abs(d))

    def make_rec(self, V):
        ix = self.ix
        Wl = self.Wl
        Wr = self.Wr
        eps = Const.EPS
        alpha = 4
        b1 = 1 / 60
        b2 = 4 / 3
        L = self.char.L(V)
        R = self.char.R(V)
        for m in range(8):
            for i in range(4, ix+4):
                # v = sum([L[i][m][n]*V[n][i] for n in range(8)])
                # vp1 = sum([L[i][m][n]*V[n][i+1] for n in range(8)])
                # vp2 = sum([L[i][m][n]*V[n][i+2] for n in range(8)])
                # vm1 = sum([L[i][m][n]*V[n][i-1] for n in range(8)])
                # vm2 = sum([L[i][m][n]*V[n][i-2] for n in range(8)])
                v = V[m][i]
                vp1 = V[m][i+1]
                vp2 = V[m][i+2]
                vm1 = V[m][i-1]
                vm2 = V[m][i-2]
                vor = b1 * (2*vm2-13*vm1+47*v+27*vp1-3*vp2)
                vmp = v + self.mm(vp1-v, alpha*(v-vm1))
                if (vor-v)*(vor-vmp) < eps:
                    Wl[m][i] = vor
                else:
                    dm1 = vm2 - 2 * vm1 + v
                    d = vm1 - 2 * v + vp1
                    dp1 = v - 2 * vp1 + vp2
                    dm4ph = self.mm4(4*d-dp1, 4*dp1-d, d, dp1)
                    dm4mh = self.mm4(4*d-dm1, 4*dm1-d, d, dm1)
                    vul = v + alpha * (v-vm1)
                    vav = 0.5 * (v+vp1)
                    vmd = vav - 0.5 * dm4ph
                    vlc = v + 0.5 * (v-vm1) + b2 * dm4mh
                    vmin = max(min(v, vp1, vmd), min(v, vul, vlc))
                    vmax = min(max(v, vp1, vmd), max(v, vul, vlc))
                    Wl[m][i] = vor + self.mm(vmin-vor, vmax-vor)
                vor = b1 * (2*vp2-13*vp1+47*v+27*vm1-3*vm2)
                vmp = v + self.mm(vm1-v, alpha*(v-vp1))
                if (vor-v)*(vor-vmp) < eps:
                    Wr[m][i-1] = vor
                else:
                    dm1 = vm2 - 2 * vm1 + v
                    d = vm1 - 2 * v + vp1
                    dp1 = v - 2 * vp1 + vp2
                    dm4ph = self.mm4(4*d-dp1, 4*dp1-d, d, dp1)
                    dm4mh = self.mm4(4*d-dm1, 4*dm1-d, d, dm1)
                    vul = v + alpha * (v-vp1)
                    vav = 0.5 * (v+vm1)
                    vmd = vav + 0.5 * dm4mh
                    vlc = v + 0.5 * (v-vp1) + b2 * dm4ph
                    vmin = max(min(v, vm1, vmd), min(v, vul, vlc))
                    vmax = min(max(v, vm1, vmd), max(v, vul, vlc))
                    Wr[m][i-1] = vor + self.mm(vmin-vor, vmax-vor)
        # Vl = self.Wl
        # Vr = self.Wr
        # for m in range(8):
        #     for i in range(4, ix+4):
        #         Vl[m][i] = sum([R[i][m][n]*Wl[n][i] for n in range(8)])
        #         Vr[m][i-1] = sum([R[i][m][n]*Wr[n][i-1] for n in range(8)])
        # return Vl, Vr
        return Wl, Wr
