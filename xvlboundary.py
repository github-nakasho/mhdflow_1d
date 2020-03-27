#!/usr/bin/env/ python3

def XVlFreeBoundary(Vl, Vr, order):
    for m in range(8):
        Vl[m][order-2] = Vr[m][order-2]
    return Vl
