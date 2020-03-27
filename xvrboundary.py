#!/usr/bin/env/ python3

def XVrFreeBoundary(Vl, Vr, ix, order):
    for m in range(8):
        Vr[m][ix+(order-2)] = Vl[m][ix+(order-2)]
    return Vr
