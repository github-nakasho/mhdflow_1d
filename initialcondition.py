#!/usr/bin/env/ python3

import numpy as np
from saveinitial import SaveInitial
from vtou import VtoU
from xulboundary import XUlFreeBoundary
from xurboundary import XUrFreeBoundary


class InitialCondition:
    def __init__(self, ixmax, order):
        self.V = np.zeros((8, ixmax))
        self.ixmax = ixmax
        self.order = order

    def RJ2a(self):
        V = self.V
        ixmax = self.ixmax
        order = self.order
        for i in range(int(ixmax/2)):
            V[0][i]=1.08
            V[1][i]=1.2
            V[2][i]=0.01
            V[3][i]=0.5
            V[4][i]=2.0/np.sqrt(4*np.pi)
            V[5][i]=3.6/np.sqrt(4*np.pi)
            V[6][i]=2.0/np.sqrt(4*np.pi)
            V[7][i]=0.95
        for i in range(int(ixmax/2), ixmax):
            V[0][i]=1.0
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=2.0/np.sqrt(4*np.pi)
            V[5][i]=4.0/np.sqrt(4*np.pi)
            V[6][i]=2.0/np.sqrt(4*np.pi)
            V[7][i]=1.0
        U = VtoU(V)
        U = XUlFreeBoundary(U, order)
        U = XUrFreeBoundary(U, ixmax, order)
        Uinitial, Vinitial = SaveInitial(U, V)
        return U, V

    def BW(self):
        V = self.V
        ixmax = self.ixmax
        order = self.order
        for i in range(int(ixmax/2)):
            V[0][i]=1.0
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=1.0
            V[4][i]=0.75
            V[5][i]=1.0
            V[6][i]=0.0
            V[7][i]=1.0
        for i in range(int(ixmax/2), ixmax):
            V[0][i]=0.125
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=-1.0
            V[4][i]=0.75
            V[5][i]=0.0
            V[6][i]=0.0
            V[7][i]=0.1
        U = VtoU(V)
        U = XUlFreeBoundary(U, order)
        U = XUrFreeBoundary(U, ixmax, order)
        Uinitial, Vinitial = SaveInitial(U, V)
        return U, V
