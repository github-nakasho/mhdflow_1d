#!/usr/bin/env/ python3

def XLeftFreeBoundary(U, order):
    for m in range(8):
        for i in range(order-1, 0, -1):
            U[m][i-1] = U[m][i]
    return U
