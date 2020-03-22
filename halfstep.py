#!/usr/bin/env/ python3

from hlld import HLLD
from reconstruct import Minmod
from source import Source
from utov import UtoV
from xleftboundary import XLeftFreeBoundary
from xrightboundary import XRightFreeBoundary

def HalfStep(U, V, dt, ixmax, order):
    Vl, Vr = Minmod()
    F = HLLD()
    S = Source()
    for m in range(8):
		for i in range(order-1, ixmax-(order-1))
            U1[m][i] = U[m] - dtodx * (F[m][i]-F[m][i-1]) + dt * S[m][i]
    U1 = XLeftFreeBoundary(U1, order)
    U1 = XRightFreeBoundary(U1, ixmax, order)
    V = UtoV(U1, V)
    return U1, V