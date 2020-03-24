#!/usr/bin/env/ python3

import numpy as np

from uvtof import UVtoF
from vtou import VtoU

def HLL(U, V, Vl, Vr, ixmax, order):
    F = np.zeros(U.shape)
    gm = 5.0 / 3.0
    # left side ----------
    vl2 = Vl[1] * Vl[1] + Vl[2] * Vl[2] + Vl[3] * Vl[3]
    bl2 = Vl[4] * Vl[4] + Vl[5] * Vl[5] + Vl[6] * Vl[6]
    vlbl = Vl[1] * Vl[4] + Vl[2] * Vl[5] + Vl[3] * Vl[6]
    Ul = VtoU(Vl)
    Fl = UVtoF(Ul, Vl, bl2, vlbl)
    gmprl = gm * Vl[7]
    vfl = np.sqrt(((bl2+gmprl)+np.sqrt((bl2+gmprl)*(bl2+gmprl)-4.0*gmprl*Vl[4]*Vl[4]))/(2.0*Vl[0]))
    # ---------- left side
    # right side ----------
    vr2 = Vr[1] * Vr[1] + Vr[2] * Vr[2] + Vr[3] * Vr[3]
    br2 = Vr[4] * Vr[4] + Vr[5] * Vr[5] + Vr[6] * Vr[6]
    vrbr = Vr[1] * Vr[4] + Vr[2] * Vr[5] + Vr[3] * Vr[6]
    Ur = VtoU(Vr)
    Fr = UVtoF(Ur, Vr, br2, vrbr)
    gmprr = gm * Vr[7]
    vfr = np.sqrt(((br2+gmprr)+np.sqrt((br2+gmprr)*(br2+gmprr)-4.0*gmprr*Vr[4]*Vr[4]))/(2.0*Vr[0]))
    # ----------- right side
    # propagation speed of Riemann fan ----------
    sl = np.minimum(Vl[1], Vr[1]) - np.maximum(vfl, vfr)
    sr = np.minimum(Vl[1], Vr[1]) + np.maximum(vfl, vfr)
    # ---------- propagation speed of Riemann fan
    # compute HLL flux
    for i in range(ixmax-(order-1)):
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
