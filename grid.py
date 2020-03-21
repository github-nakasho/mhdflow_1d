#!/usr/bin/env/ python3

import numpy as np

def Grid(ix, order, x1, x0):
    ixmax = ix + 2 * (order-1)
    dx = (x1-x0) / (ix-1)
    x = np.zeros(ixmax)
    x[0] = - dx
    for i in range(ixmax-1):
        x[i+1] = x[i] + dx
    minlength = dx
    return {'ixmax': ixmax, 
            'minlength': minlength, 
            'dx': dx, 
            'x': x}
    