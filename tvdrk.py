#!/usr/bin/env/ python3

from fullstep import FullStep
from halfstep import HalfStep

def TVDRK(U, V, ixmax, dx, dt, order):
    U1, V = HalfStep(U, V, ixmax, dx, dt, order)
    U, V = FullStep(U, U1, V, ixmax, dx, dt, order)
    return U, V
