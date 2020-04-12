#!/usr/bin/env/ python3

import numpy as np

from convert import Convert


class InitialCondition:
    def __init__(self, x, ix, order):
        self.x = x
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
            V[4][i]=0.95
            V[5][i]=2.0/np.sqrt(4*np.pi)
            V[6][i]=3.6/np.sqrt(4*np.pi)
            V[7][i]=2.0/np.sqrt(4*np.pi)
        for i in range(int(ixmax/2), ixmax):
            V[0][i]=1.0
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=1.0
            V[5][i]=2.0/np.sqrt(4*np.pi)
            V[6][i]=4.0/np.sqrt(4*np.pi)
            V[7][i]=2.0/np.sqrt(4*np.pi)
        U = self.conv.VtoU(V)
        return U, V

    def BW(self):
        V = self.V
        ixmax = self.ixmax
        for i in range(int(ixmax/2)):
            V[0][i]=1.0
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=1.0
            V[5][i]=0.75
            V[6][i]=1.0
            V[7][i]=0.0
        for i in range(int(ixmax/2), ixmax):
            V[0][i]=0.125
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=0.1
            V[5][i]=0.75
            V[6][i]=-1.0
            V[7][i]=0.0
        U = self.conv.VtoU(V)
        return U, V

    def SO(self):
        V = self.V
        ixmax = self.ixmax
        for i in range(int(ixmax/10)):
            V[0][i]=3.857143
            V[1][i]=2.629369
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=0.0
            V[5][i]=0.0
            V[6][i]=0.0
            V[7][i]=10.33333
        for i in range(int(ixmax/10), ixmax):
            V[0][i]=1 + 0.2 * np.sin(5*np.pi*self.x[i])
            V[1][i]=0.0
            V[2][i]=0.0
            V[3][i]=0.0
            V[4][i]=0.0
            V[5][i]=0.0
            V[6][i]=0.0
            V[7][i]=1.0
        U = self.conv.VtoU(V)
        return U, V
