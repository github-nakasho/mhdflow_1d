#!/usr/bin/env/ python3

def XVrFreeBoundary(Vl, Vr, ixmax, order):
    for m in range(8):
        Vr[m][ixmax-(order-1)-1] = Vl[m][ixmax-(order-1)-1]
    return Vr
