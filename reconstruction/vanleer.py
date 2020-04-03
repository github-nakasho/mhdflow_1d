#!/usr/bin/env/ python3

import numpy as np


class VanLeer:
    def __init__(self, ix, order):
        ixmax = ix + 2 * (order-1)
        self.ix = ix
        self.Wl = np.zeros((8, ixmax-1))
        self.Wr = np.zeros((8, ixmax-1))
        
    def make_rec(self, V):
        ix = self.ix
        Wl = self.Wl
        Wr = self.Wr
        for m in range(8):
            for i in range(1, ix+1):
                w = V[m][i]
                wp = V[m][i+1]
                wm = V[m][i-1]
                a = w - wm
                b = wp - w
                if a * b <= 0:
                    grad = 0
                else:
                    grad = a * b / (a+b)
                Wl[m][i] = w + grad
                Wr[m][i-1] = w - grad
        return Wl, Wr
