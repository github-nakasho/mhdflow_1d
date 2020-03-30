#!/usr/bin/env/ python3

import numpy as np

from convert import Convert


class InitialCondition:
    def __init__(self, ix, order):
        self.ixmax = ix + 2 * (order-1)
        self.V = np.zeros((8, self.ixmax))
        self.conv = Convert()

    def RJ2a(self):
        V = self.V
        ixmax = self.ixmax
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
        U = self.conv.VtoU(V)
        return U, V

    def BW(self):
        V = self.V
        ixmax = self.ixmax
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
        U = self.conv.VtoU(V)
        return U, V
