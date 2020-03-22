#!/usr/bin/env/ python3

from fullstep import FullStep
from halfstep import HalfStep

def TVDRK(U, V, dt):
    U, V = HalfStep(U, V, dt)
    U, V = FullStep(U, V, dt)
    return U, V
