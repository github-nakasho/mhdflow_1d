#!/usr/bin/env/ python3

def XVlFreeBoundary(Vl, Vr, order):
    for m in range(8):
        Vl[m][0] = Vr[m][1]
    return Vl
