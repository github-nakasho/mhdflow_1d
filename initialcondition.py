#!/usr/bin/env/ python3

import numpy as np

def DW1aShockTube(ixmax):
    V = np.zeros((8, ixmax))
    U = np.zeros((8, ixmax))
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
#     U = VtoU()
#     XLeftFreeBoundary(U)
#     XRightFreeBoundary(U)
#     SaveInitial()
    return U, V
