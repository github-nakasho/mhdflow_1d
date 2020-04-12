#!/usr/bin/env/ python3

import numpy as np

from common import Const


class Characteristic:
    def __init__(self, ix, order):
        self.ix = ix
        self.order = order

    def make_coeff(self, V, i):
        gm = Const.GAMMA
        ro = V[0][i]
        roi = 1 / ro
        sro = np.sqrt(ro)
        gmpr = gm * V[4][i]
        a = np.sqrt(gmpr*roi)
        bx = V[5][i]
        bx2 = bx ** 2
        bt2 = V[6][i] ** 2 + V[7][i] ** 2
        ca = np.sqrt((bx2+bt2)*roi)
        cax = np.sqrt(bx2*roi)
        aca = a * a + ca * ca
        cf = np.sqrt(0.5*(aca+np.sqrt(aca**2-4*a*a*cax*cax)))
        cs = np.sqrt(0.5*(aca-np.sqrt(aca**2-4*a*a*cax*cax)))
        if cf * cf - cs * cs <= 0:
            alf = 1
            als = 0
        elif a * a - cs * cs <= 0:
            alf = 0
            als = 1
        elif cf * cf - a * a <= 0:
            alf = 1.0
            als = 0.0
        else:
            alf = np.sqrt((a*a-cs*cs)/(cf*cf-cs*cs))
            als = np.sqrt((cf*cf-a*a)/(cf*cf-cs*cs))
        s = np.sign(bx)
        nf = 0.5 / a ** 2
        ns = 0.5 / a ** 2
        cff = cf * alf
        css = cs * als
        qf = cff * s
        qs = css * s
        aaf = a * alf * sro
        aas = a * als * sro
        by = 0
        bz = 0
        if bt2 != 0:
            by = V[6][i] / np.sqrt(bt2)
            bz = V[7][i] / np.sqrt(bt2)
        return nf, cff, qs, by, bz, alf, roi, aas, s, sro, css, qf, als, aaf, a, ns, ro

    def L(self, V):
        L = [0] * len(V[0])
        for i in range(len(V[0])):
            nf, cff, qs, by, bz, alf, roi, aas, s, sro, css, qf, als, aaf, a, ns, ro = self.make_coeff(V, i)
            L[i] = np.array([nf*np.array([0, -cff, qs*by, qs*bz, alf*roi, aas*by*roi, aas*bz*roi]), 
                            0.5*np.array([0, 0, -bz, by, 0, -bz*s/sro, by*s/sro]), 
                            ns*np.array([0, -css, -qf*by, -qf*bz, als*roi, -aaf*by*roi, -aaf*bz*roi]), 
                            np.array([1, 0, 0, 0, -1/a**2, 0, 0]), 
                            ns*np.array([0, css, qf*by, qf*bz, als*roi, -aaf*by*roi, -aaf*bz*roi]), 
                            0.5*np.array([0, 0, bz, -by, 0, -bz*s/sro, by*s/sro]),
                            nf*np.array([0, cff, -qs*by, -qs*bz, alf*roi, aas*by*roi, aas*bz*roi])])
        return np.array(L)
    
    def R(self, V):
        R = [0] * len(V[0])
        for i in range(len(V[0])):
            nf, cff, qs, by, bz, alf, roi, aas, s, sro, css, qf, als, aaf, a, ns, ro = self.make_coeff(V, i)
            R[i] = np.array([np.array([ro*alf, 0, ro*als, 1, ro*als, 0, ro*alf]), 
                            np.array([-cff, 0, -css, 0, css, 0, cff]), 
                            np.array([qs*by, -bz, -qf*by, 0, qf*by, bz, -qs*by]), 
                            np.array([qs*bz, by, -qf*bz, 0, qf*bz, -by, -qs*bz]), 
                            ro*a**2*np.array([alf, 0, als, 0, als, 0, alf]), 
                            np.array([aas*by, -bz*s*sro, -aaf*by, 0, -aaf*by, -bz*s*sro, aas*by]), 
                            np.array([aas*bz, by*s*sro, -aaf*bz, 0, -aaf*bz, by*s*sro, aas*bz])])
        return np.array(R)
