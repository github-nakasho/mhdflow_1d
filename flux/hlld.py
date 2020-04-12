#!/usr/bin/env/ python3

import numpy as np

from common import Const
from convert import Convert

class HLLD:
    def __init__(self):
        self.conv = Convert()

    def make_flux(self, Vl, Vr, ix, order):
        eps = Const.EPS
        gm = Const.GAMMA
        # set primitive variables @ left-face
        rol = Vl[0]
        vnl = Vl[1]
        vtl = Vl[2]
        vul = Vl[3]
        prl = Vl[4]
        bnc = Vl[5]
        btl = Vl[6]
        bul = Vl[7]
        # set primitive variables @ right-face
        ror = Vr[0]
        vnr = Vr[1]
        vtr = Vr[2]
        vur = Vr[3]
        prr = Vr[4]
        bnc = Vr[5]
        btr = Vr[6]
        bur = Vr[7]
        # bn @ the interface
        bnc2 = bnc * bnc
        sgn = np.sign(bnc)
        # variables @ the left-face
        roli = 1 / rol
        pml = 0.5 * (btl**2+bul**2)
        ptl = prl + pml
        enl = gm * prl + pml + 0.5 * rol * (vnl**2+vtl**2+vul**2)
        vbl = vtl * btl + vul * bul
        # variables @ the right-face
        rori = 1 / ror
        pmr = 0.5 * (btr**2+bur**2)
        ptr = prr + pmr
        enr = gm * prr + pmr + 0.5 * ror * (vnr**2+vtr**2+vur**2)
        vbr = vtr * btr + vur * bur
        # maximum / minimum wave speeds
        cl2 = gm * prl * roli
        cr2 = gm * prr * rori
        cal2 = bnc2 * roli
        car2 = bnc2 * rori
        cbl2 = cl2 + cal2 + 2 * pml * roli
        cbr2 = cr2 + car2 + 2 * pmr * rori
        cfl2 = 0.5 * (cbl2+np.sqrt(abs(cbl2*cbl2-4*cl2*cal2)))
        cfr2 = 0.5 * (cbr2+np.sqrt(abs(cbr2*cbr2-4*cr2*car2)))
        cfl = np.sqrt(cfl2)
        cfr = np.sqrt(cfr2)
        sl = np.minimum(0.0, np.minimum(vnl, vnr)-np.maximum(cfl, cfr))
        sr = np.maximum(0.0, np.maximum(vnl, vnr)+np.maximum(cfl, cfr))
        # HLL average of the normal velocity and the total pressure
        slvl = sl - vnl
        srvr = sr - vnr
        rslvl = rol * slvl
        rsrvr = ror * srvr
        drsvi = 1 / (rsrvr-rslvl)
        vnc = (rsrvr*vnr-rslvl*vnl-ptr+ptl) * drsvi
        ptc = (rsrvr*ptl-rslvl*ptr+rsrvr*rslvl*(vnr-vnl)) * drsvi
        # variables of the outer sides in the Riemann fan
        slvc = sl - vnc
        srvc = sr - vnc
        ro2l = rslvl / slvc
        ro2r = rsrvr / srvc
        rhdl = rslvl * slvc - bnc2
        rhdr = rsrvr * srvc - bnc2
        # set ta (= true array) & fa (= false array)
        ta = abs(rhdl) > eps
        rhdli = np.zeros(np.shape(rhdl))
        rhnvl = np.zeros(np.shape(rhdl))
        rhnbl = np.zeros(np.shape(rhdl))
        vt2l = np.zeros(np.shape(rhdl))
        vu2l = np.zeros(np.shape(rhdl))
        bt2l = np.zeros(np.shape(rhdl))
        bu2l = np.zeros(np.shape(rhdl))
        # if abs(rhdl) > eps
        rhdli[ta] = 1 / rhdl[ta]
        rhnvl[ta] = (vnl[ta]-vnc[ta]) * bnc[ta]
        rhnbl[ta] = rslvl[ta] * slvl[ta] - bnc2[ta]
        vt2l[ta] = vtl[ta] + rhnvl[ta] * rhdli[ta] * btl[ta]
        vu2l[ta] = vul[ta] + rhnvl[ta] * rhdli[ta] * bul[ta]
        bt2l[ta] = rhnbl[ta] * rhdli[ta] * btl[ta]
        bu2l[ta] = rhnbl[ta] * rhdli[ta] * bul[ta]
        # else 
        fa = ta == False
        vt2l[fa] = vtl[fa]
        vu2l[fa] = vul[fa]
        bt2l[fa] = btl[fa]
        bu2l[fa] = bul[fa]
        # set ta (= true array) & fa (= false array)
        ta = abs(rhdr) > eps
        rhdri = np.zeros(np.shape(rhdr))
        rhnvr = np.zeros(np.shape(rhdr))
        rhnbr = np.zeros(np.shape(rhdr))
        vt2r = np.zeros(np.shape(rhdr))
        vu2r = np.zeros(np.shape(rhdr))
        bt2r = np.zeros(np.shape(rhdr))
        bu2r = np.zeros(np.shape(rhdr))
        # if abs(rhdr) > eps
        rhdri[ta] = 1 / rhdr[ta]
        rhnvr[ta] = (vnr[ta]-vnc[ta]) * bnc[ta]
        rhnbr[ta] = rsrvr[ta] * srvr[ta] - bnc2[ta]
        vt2r[ta] = vtr[ta] + rhnvr[ta] * rhdri[ta] * btr[ta]
        vu2r[ta] = vur[ta] + rhnvr[ta] * rhdri[ta] * bur[ta]
        bt2r[ta] = rhnbr[ta] * rhdri[ta] * btr[ta]
        bu2r[ta] = rhnbr[ta] * rhdri[ta] * bur[ta]
        # else
        fa = ta & False
        vt2r[fa] = vtr[fa]
        vu2r[fa] = vur[fa]
        bt2r[fa] = btr[fa]
        bu2r[fa] = bur[fa]
        vb2l = vt2l * bt2l + vu2l * bu2l
        vb2r = vt2r * bt2r + vu2r * bu2r
        en2l = (slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l)) / slvc
        en2r = (srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r)) / srvc
        # variables of the inner sides in Riemann fan
        rro2l = np.sqrt(ro2l)
        rro2r = np.sqrt(ro2r)
        rro2i = 1 / (rro2r+rro2l)
        vt3m = (rro2r*vt2r+rro2l*vt2l+(bt2r-bt2l)*sgn) * rro2i
        vu3m = (rro2r*vu2r+rro2l*vu2l+(bu2r-bu2l)*sgn) * rro2i
        bt3m = (rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn) * rro2i
        bu3m = (rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn) * rro2i
        vb3m = vt3m * bt3m + vu3m * bu3m
        en3l = en2l - rro2l * (vb2l-vb3m) * sgn
        en3r = en2r + rro2r * (vb2r-vb3m) * sgn
        # variables @ the interface
        ta1 = vnc - abs(bnc) / rro2l > 0
        rou = np.zeros(np.shape(ro2l))
        vtu = np.zeros(np.shape(ro2l))
        vuu = np.zeros(np.shape(ro2l))
        btu = np.zeros(np.shape(ro2l))
        buu = np.zeros(np.shape(ro2l))
        enu = np.zeros(np.shape(ro2l))
        rou[ta1] = ro2l[ta1]
        vtu[ta1] = vt2l[ta1]
        vuu[ta1] = vu2l[ta1]
        btu[ta1] = bt2l[ta1]
        buu[ta1] = bu2l[ta1]
        enu[ta1] = en2l[ta1]
        ta2 =  (vnc - abs(bnc) / rro2l <= 0) & (vnc >= 0)
        rou[ta2] = ro2l[ta2]
        vtu[ta2] = vt3m[ta2]
        vuu[ta2] = vu3m[ta2]
        btu[ta2] = bt3m[ta2]
        buu[ta2] = bu3m[ta2]
        enu[ta2] = en3l[ta2]
        ta3 = (vnc < 0.0) & (vnc + abs(bnc) / rro2r >= 0)
        rou[ta3] = ro2r[ta3]
        vtu[ta3] = vt3m[ta3]
        vuu[ta3] = vu3m[ta3]
        btu[ta3] = bt3m[ta3]
        buu[ta3] = bu3m[ta3]
        enu[ta3] = en3r[ta3]
        ta4 = ta1 | ta2 | ta3
        fa = ta4 & False
        rou[fa] = ro2r[fa]
        vtu[fa] = vt2r[fa]
        vuu[fa] = vu2r[fa]
        btu[fa] = bt2r[fa]
        buu[fa] = bu2r[fa]
        enu[fa] = en2r[fa]
        # HLLD fluxes
        F = self.conv.PtoF(rou, vnc, vtu, vuu, bnc, btu, buu, ptc, enu)
        return F
