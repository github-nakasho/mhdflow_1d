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

def FullStep(U, U1, V, ixmax, dx, dt, order):
    Vl, Vr = Minmod(V, ixmax)
    Vl = XVlFreeBoundary(Vl, Vr, order)
    Vr = XVrFreeBoundary(Vl, Vr, ixmax, order)
    F = HLL(U1, V, Vl, Vr, ixmax, order)
    S = Source(U1)
    for m in range(8):
        for i in range(order-1, ixmax-(order-1)):
            U[m][i] = 0.5 * (U[m][i]+U1[m][i]-dt/dx*(F[m][i]-F[m][i-1])+dt*S[m][i])
    U = XUlFreeBoundary(U, order)
    U = XUrFreeBoundary(U, ixmax, order)
    V = UtoV(U, V)
    return U, V