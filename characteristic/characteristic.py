#!/usr/bin/env/ python3

import numpy as np

from common import Const


class Characteristic:
    def __init__(self, ix, order):
        self.ix = ix
        self.order = order

    def L(self, V):
        gm = Const.GAMMA
        L = [0] * len(V[0])
        for i in range(len(V[0])):
            roi = 1 / V[0][i]
            sro = np.sqrt(V[0][i])
            a = np.sqrt(gm*V[7][i]*roi)
            ca = np.sqrt((V[4][i]**2+V[5][i]**2+V[6][i]**2)*roi)
            cax = np.sqrt(V[4][i]**2*roi)
            cf = np.sqrt(0.5*(a*a+ca*ca+np.sqrt((a*a+ca*ca)**2-4*a*a*cax*cax)))
            cs = np.sqrt(0.5*(a*a+ca*ca-np.sqrt((a*a+ca*ca)**2-4*a*a*cax*cax)))
            alf = np.sqrt((a*a-cs*cs)/(cf*cf-cs*cs))
            als = np.sqrt((cf*cf-a*a)/(cf*cf-cs*cs))
            cff = cf * alf
            css = cs * als
            s = np.sign(V[4][i])
            qf = cf * alf * s
            qs = cs * als * s
            aaf = a * alf * sro
            aas = a * als * sro
            by = V[5][i] / np.sqrt(V[5][i]**2+V[6][i]**2)
            bz = V[6][i] / np.sqrt(V[5][i]**2+V[6][i]**2)
            nf = 1 / (2*a*a)
            ns = 1 / (2*a*a)
            L[i] = np.array([[0, -nf*cff, nf*qs*by, nf*qs*bz, 0, nf*aas*by*roi, nf*aas*bz*roi, nf*alf*roi], 
                            [0, 0, -bz/2, by/2, 0, -bz*s/(2*sro), by*s/(2*sro), 0], 
                            [0, -ns*css, -ns*qf*by, -ns*qf*bz, 0, -ns*aaf*by*roi, -ns*aaf*bz*roi, ns*als*roi], 
                            [1, 0, 0, 0, 0, 0, 0, -1/(a*a)], 
                            [0, 0, 0, 0, 1, 0, 0, 0], 
                            [0, 0, bz/2, -by/2, 0, -bz*s/(2*sro), by*s/(2*sro), 0],
                            [0, nf*cff, -nf*qs*by, -nf*qs*bz, 0, nf*aas*by*roi, nf*aas*bz*roi, nf*alf*roi], 
                            [0, ns*css, ns*qf*by, ns*qf*bz, 0, -ns*aaf*by*roi, -ns*aaf*bz*roi, ns*als*roi]])
        return np.array(L)
    
    def R(self, V):
        gm = Const.GAMMA
        R = [0] * len(V[0])
        for i in range(len(V[0])):
            ro = V[0][i]
            roi = 1 / ro
            sro = np.sqrt(V[0][i])
            a = np.sqrt(gm*V[7][i]*roi)
            ca = np.sqrt((V[4][i]**2+V[5][i]**2+V[6][i]**2)*roi)
            cax = np.sqrt(V[4][i]**2*roi)
            cf = np.sqrt(0.5*(a*a+ca*ca+np.sqrt((a*a+ca*ca)**2-4*a*a*cax*cax)))
            cs = np.sqrt(0.5*(a*a+ca*ca-np.sqrt((a*a+ca*ca)**2-4*a*a*cax*cax)))
            alf = np.sqrt((a*a-cs*cs)/(cf*cf-cs*cs))
            als = np.sqrt((cf*cf-a*a)/(cf*cf-cs*cs))
            cff = cf * alf
            css = cs * als
            s = np.sign(V[4][i])
            qf = cf * alf * s
            qs = cs * als * s
            aaf = a * alf * sro
            aas = a * als * sro
            by = V[5][i] / np.sqrt(V[5][i]**2+V[6][i]**2)
            bz = V[6][i] / np.sqrt(V[5][i]**2+V[6][i]**2)
            R[i] = np.array([[ro*alf, 0, ro*als, 1, 0, 0, ro*alf, ro*als], 
                            [-cff, 0, -css, 0, 0, 0, cff, css], 
                            [qs*by, -bz, -qf*by, 0, 0, bz, -qs*by, qf*by], 
                            [qs*bz, by, -qf*bz, 0, 0, -by, -qs*bz, qf*bz], 
                            [0, 0, 0, 0, 1, 0, 0, 0], 
                            [aas*by, -bz*s*sro, -aaf*by, 0, 0, -bz*s*sro, aas*by, -aaf*by], 
                            [aas*bz, by*s*sro, -aaf*bz, 0, 0, by*s*sro, aas*bz, -aaf*bz], 
                            [ro*a*a*alf, 0, ro*a*a*als, 0, 0, 0, ro*a*a*alf, ro*a*a*als]])
        return np.array(R)
