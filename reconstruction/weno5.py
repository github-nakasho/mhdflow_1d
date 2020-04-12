#!/usr/bin/env/ python3

import numpy as np

from characteristic import Characteristic
from common import Const


class WENO5:
    def __init__(self, ix, order):
        order = 5
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Vl = np.zeros((8, ixmax-1))
        self.Vr = np.zeros((8, ixmax-1))
        

    def make_rec(self, V):
        ix = self.ix
        Vl = self.Vl
        Vr = self.Vr
        eps = Const.EPS
        cd = np.array([1/10, 6/10, 3/10])
        c0 = np.array([2/6, -7/6, 11/6])
        c1 = np.array([-1/6, 5/6, 2/6])
        c2 = np.array([2/6, 5/6, -1/6])
        isc1 = 13 / 12
        isc2 = 1 / 4
        iis = np.zeros(3)
        a = np.zeros(3)
        om = np.zeros(3)
        for m in range(8):
            for i in range(4, ix+4):
                v = V[m][i]
                vp1 = V[m][i+1]
                vp2 = V[m][i+2]
                vm1 = V[m][i-1]
                vm2 = V[m][i-2]
                iis[0] = isc1 * (vm2-2*vm1+v) ** 2 + isc2 * (vm2-4*vm1+3*v) ** 2
                iis[1] = isc1 * (vm1-2*v+vp1) ** 2 + isc2 * (vm1-vp1) ** 2
                iis[2] = isc1 * (v-2*vp1+vp2) ** 2 + isc2 * (3*v-4*vp1+vp2) ** 2 
                a = cd / (iis+eps) ** 2
                om = a / sum(a)
                v0 = c0[0] * vm2 + c0[1] * vm1 + c0[2] * v
                v1 = c1[0] * vm1 + c1[1] * v + c1[2] * vp1
                v2 = c2[0] * v + c2[1] * vp1 + c2[2] * vp2
                Vl[m][i] = om[0] * v0 + om[1] * v1 + om[2] * v2
                iis[0] = isc1 * (vp2-2*vp1+v) ** 2 + isc2 * (vp2-4*vp1+3*v) ** 2
                iis[1] = isc1 * (vp1-2*v+vm1) ** 2 + isc2 * (vp1-vm1) ** 2
                iis[2] = isc1 * (v-2*vm1+vm2) ** 2 + isc2 * (3*v-4*vm1+vm2) ** 2 
                a = cd / (iis+eps) ** 2
                om = a / sum(a)
                v0 = c0[0] * vp2 + c0[1] * vp1 + c0[2] * v
                v1 = c1[0] * vp1 + c1[1] * v + c1[2] * vm1
                v2 = c2[0] * v + c2[1] * vm1 + c2[2] * vm2
                Vr[m][i-1] = om[0] * v0 + om[1] * v1 + om[2] * v2
        return Vl, Vr
        