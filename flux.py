#!/usr/bin/env/ python3

import numpy as np

from uvtof import UVtoF
from vtou import VtoU


class Flux:
    def __init__(self):
        pass

    def HLL(self, Vl, Vr, ix, order):
        F = np.zeros(Vl.shape)
        gm = 5.0 / 3.0
        # left side ----------
        bl2 = Vl[4] * Vl[4] + Vl[5] * Vl[5] + Vl[6] * Vl[6]
        vlbl = Vl[1] * Vl[4] + Vl[2] * Vl[5] + Vl[3] * Vl[6]
        Ul = VtoU(Vl)
        Fl = UVtoF(Ul, Vl, bl2, vlbl)
        gmprl = gm * Vl[7]
        vfl = np.sqrt(((bl2+gmprl)+np.sqrt((bl2+gmprl)*(bl2+gmprl)-4*gmprl*Vl[4]*Vl[4]))/(2*Vl[0]))
        # ---------- left side
        # right side ----------
        br2 = Vr[4] * Vr[4] + Vr[5] * Vr[5] + Vr[6] * Vr[6]
        vrbr = Vr[1] * Vr[4] + Vr[2] * Vr[5] + Vr[3] * Vr[6]
        Ur = VtoU(Vr)
        Fr = UVtoF(Ur, Vr, br2, vrbr)
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

    def HLLD(self, Vl, Vr, ix, order):
        gm = 5 / 3
        vl2 = Vl[1] * Vl[1] + Vl[2] * Vl[2] + Vl[3] * Vl[3]
        bl2 = Vl[4] * Vl[4] + Vl[5] * Vl[5] + Vl[6] * Vl[6]
        vlbl = Vl[1] * Vl[4] + Vl[2] * Vl[5] + Vl[3] * Vl[6]
        ptl = Vl[7] + 0.5 * bl2
        gmprl = gm * Vl[7]
        vfl = np.sqrt(((bl2+gmprl)+np.sqrt((bl2+gmprl)*(bl2+gmprl)-4*gmprl*Vl[4]*Vl[4]))/(2*Vl[0]))
        vr2 = Vr[1] * Vr[1] + Vr[2] * Vr[2] + Vr[3] * Vr[3]
        br2 = Vr[4] * Vr[4] + Vr[5] * Vr[5] + Vr[6] * Vr[6]
        vrbr = Vr[1] * Vr[4] + Vr[2] * Vr[5] + Vr[3] * Vr[6]
        ptr = Vr[7] + 0.5 * br2
        gmprr = gm * Vr[7]
        vfr = np.sqrt(((br2+gmprr)+np.sqrt((br2+gmprr)*(br2+gmprr)-4*gmprr*Vr[4]*Vr[4]))/(2*Vr[0]))
        # Riemann fan speed ------
        sl = np.minimum(Vl[1], Vr[1]) - np.maximum(vfl, vfr)
        sr = np.maximum(Vl[1], Vr[1]) + np.maximum(vfl, vfr)
        if sl > 0.0:
            Ul = VtoU(Vl)
            Fl = UVtoF(Ul, Vl, bl2, vlbl)
        elif sr < 0.0:
            Ur = VtoU(Vr)
            Fr = UVtoF(Ur, Vr, br2, vrbr)
        else:
            Ul = VtoU(Vl)
            Ur = VtoU(Vr)
            slvl = sl - Vl[1]
            srvr = sr - Vr[1]
            slvlrhol = slvl * Ul[0]
            srvrrhor = srvr * Ur[0]
            sm = (srvr*Ur[1]-slvl*Ul[1]-ptr+ptl)/(srvrrhor-slvlrhol)
            slsm = sl - sm
            srsm = sr - sm
            oneoslsm = 1 / slsm
            oneosrsm = 1 / srsm
            ptint=(srvrrhor*ptl-slvlrhol*ptr+srvrrhor*slvlrhol*(Vr[1]-Vl[1])) / (srvrrhor-slvlrhol)
            bxint = (sr*Ur[4]-sl*Ul[4]) / (sr-sl)
            bxint2 = bxint * bxint
            rholint = slvlrhol * oneoslsm
            rhorint = srvrrhor * oneosrsm
            sign = abs(abs(slvlrhol*slsm-bxint2)-epsilon*bxint2) / (abs(slvlrhol*slsm-bxint2)-epsilon*bxint2)
            den = ((1.0/(slvlrhol*slsm-bxint2)+1.0) + (1.0/(slvlrhol*slsm-bxint2)-1.0)*sign) * 0.5
            num = (((slvlrhol*slvl-bxint2)+1.0) + ((slvlrhol*slvl-bxint2)-1.0)*sign) * 0.5
            vylint = Vl[2] - (bxint*Ul[5]*(sm-Vl[1])*(den+den*sign)*0.5)
            vzlint = Vl[3] - (bxint*Ul[6]*(sm-Vl[1])*(den+den*sign)*0.5)
            bylint = Ul[5] * num * den
            bzlint = Ul[6] * num * den
            vlintblint = sm * bxint + vylint * bylint + vzlint * bzlint
            sign = abs(abs(srvrrhor*srsm-bxint2)-epsilon*bxint2) / (abs(srvrrhor*srsm-bxint2)-epsilon*bxint2)
            den = ((1.0/(srvrrhor*srsm-bxint2)+1.0) + (1.0/(srvrrhor*srsm-bxint2)-1.0)*sign) * 0.5
            num = (((srvrrhor*srvr-bxint2)+1.0) + ((srvrrhor*srvr-bxint2)-1.0)*sign) * 0.5
            vyrint = Vr[2] - (bxint*Ur[5]*(sm-Vr[1])*(den+den*sign)*0.5)
            vzrint = Vr[3] - (bxint*Ur[6]*(sm-Vr[1])*(den+den*sign)*0.5)
            byrint = Ur[5] * num * den
            bzrint = Ur[6] * num * den
            vrintbrint = sm * bxint + vyrint * byrint + vzrint * bzrint
            elint = (slvl*Ul[7]-ptl*Vl[1]+ptint*sm+bxint*(vlbl-vlintblint))*oneoslsm
            erint = (srvr*Ur[7]-ptr*Vr[1]+ptint*sm+bxint*(vrbr-vrintbrint))*oneosrsm
            sqrtrholint = np.sqrt(rholint)
            sqrtrhorint = np.sqrt(rhorint)
            slint = sm - abs(bxint) / sqrtrholint
            srint = sm + abs(bxint) / sqrtrhorint
            if sm >= 0:
                if slint >= 0:
                    PtoF()
                else:
