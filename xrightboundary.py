#!/usr/bin/env/ python3

def XRightFreeBoundary(U, ixmax, order):
    for m in range(8):
        for i in range(ixmax-(order-1), ixmax):
            U[m][i] = U[m][i-1]
    return U
