#!/usr/bin/env/ python3

import numpy as np

from common import Const
from convert import Convert


class HLL:
    def __init__(self):
        self.convert = Convert

    def make_flux(self, Vl, Vr, ix, order):
        F = np.zeros(Vl.shape)
        gm = Const.GAMMA
        # left side ----------
        bl2 = Vl[4] * Vl[4] + Vl[5] * Vl[5] + Vl[6] * Vl[6]
        vlbl = Vl[1] * Vl[4] + Vl[2] * Vl[5] + Vl[3] * Vl[6]
        Ul = self.convert.VtoU(Vl)
        Fl = self.convert.UVtoF(Ul, Vl, bl2, vlbl)
        gmprl = gm * Vl[7]
        vfl = np.sqrt(((bl2+gmprl)+np.sqrt((bl2+gmprl)*(bl2+gmprl)-4*gmprl*Vl[4]*Vl[4]))/(2*Vl[0]))
        # ---------- left side
        # right side ----------
        br2 = Vr[4] * Vr[4] + Vr[5] * Vr[5] + Vr[6] * Vr[6]
        vrbr = Vr[1] * Vr[4] + Vr[2] * Vr[5] + Vr[3] * Vr[6]
        Ur = self.convert.VtoU(Vr)
        Fr = self.convert.UVtoF(Ur, Vr, br2, vrbr)
        gmprr = gm * Vr[7]
        vfr = np.sqrt(((br2+gmprr)+np.sqrt((br2+gmprr)*(br2+gmprr)-4*gmprr*Vr[4]*Vr[4]))/(2*Vr[0]))
        # ----------- right side
        # propagation speed of Riemann fan ----------
        sl = np.minimum(Vl[1], Vr[1]) - np.maximum(vfl, vfr)
        sr = np.maximum(Vl[1], Vr[1]) + np.maximum(vfl, vfr)
        # ---------- propagation speed of Riemann fan
        # compute HLL flux
        for i in range(order-2, ix+(order-2)+1):
            if sl[i] > 0:
                for m in range(8):
                    F[m][i] = Fl[m][i]
            elif sr[i] < 0:
                for m in range(8):
                    F[m][i] = Fr[m][i]
            else:
                for m in range(8):
                    F[m][i] = (sr[i]*Fl[m][i]-sl[i]*Fr[m][i]+sr[i]*sl[i]*(Ur[m][i]-Ul[m][i])) / (sr[i]-sl[i])
        return F
        