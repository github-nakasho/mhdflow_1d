#!/usr/bin/env/ python3

import numpy as np

from flux import HLL
from reconstruct import Minmod
from source import Source
from utov import UtoV
from xulboundary import XUlFreeBoundary
from xurboundary import XUrFreeBoundary
from xvlboundary import XVlFreeBoundary
from xvrboundary import XVrFreeBoundary

def HalfStep(U, V, ixmax, dx, dt, order):
    U1 = np.zeros(U.shape)
    Vl, Vr = Minmod(V, ixmax, order, dx, dt)
    Vl = XVlFreeBoundary(Vl, Vr, order)
    Vr = XVrFreeBoundary(Vl, Vr, ixmax, order)
    F = HLL(Vl, Vr, ixmax, order)
    S = Source(U)
    for m in range(8):
        for i in range(order-1, ixmax-(order-1)):
            U1[m][i] = U[m][i] - dt / dx * (F[m][i]-F[m][i-1]) + dt * S[m][i]
    U1 = XUlFreeBoundary(U1, order)
    U1 = XUrFreeBoundary(U1, ixmax, order)
    V = UtoV(U1)
    return U1, V