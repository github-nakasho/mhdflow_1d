#!/usr/bin/env/ python3

import numpy as np

from flux import Flux
from reconstruct import Reconstruct
from source import Source
from utov import UtoV
from xulboundary import XUlFreeBoundary
from xurboundary import XUrFreeBoundary
from xvlboundary import XVlFreeBoundary
from xvrboundary import XVrFreeBoundary

class TVDRK:
    def __init__(self, ix, dx, order):
        self.ix = ix
        self.dx = dx
        self.order = order
        self.ixmax = ix + 2 * (order-1)
        self.rec = Reconstruct(self.ix, self.order)
        self.flux = Flux()
    
    def HalfStep(self, U, V, dt):
        ix = self.ix
        dx = self.dx
        order = self.order
        ixmax = self.ixmax
        U1 = np.zeros(U.shape)
        Vl, Vr = self.rec.Minmod(V)
        Vl = XVlFreeBoundary(Vl, Vr, order)
        Vr = XVrFreeBoundary(Vl, Vr, ix, order)
        F = self.flux.HLLD(Vl, Vr, ix, order)
        S = Source(U)
        for m in range(8):
            for i in range(order-1, ix+(order-1)):
                U1[m][i] = U[m][i] - dt / dx * (F[m][i]-F[m][i-1]) + dt * S[m][i]
        U1 = XUlFreeBoundary(U1, order)
        U1 = XUrFreeBoundary(U1, ixmax, order)
        V = UtoV(U1)
        return U1, V

    def FullStep(self, U, U1, V, dt):
        ix = self.ix
        dx = self.dx
        order = self.order
        ixmax = self.ixmax
        Vl, Vr = self.rec.Minmod(V)
        Vl = XVlFreeBoundary(Vl, Vr, order)
        Vr = XVrFreeBoundary(Vl, Vr, ix, order)
        F = self.flux.HLLD(Vl, Vr, ix, order)
        S = Source(U1)
        for m in range(8):
            for i in range(order-1, ix+(order-1)):
                U[m][i] = 0.5 * (U[m][i]+U1[m][i]-dt/dx*(F[m][i]-F[m][i-1])+dt*S[m][i])
        U = XUlFreeBoundary(U, order)
        U = XUrFreeBoundary(U, ixmax, order)
        V = UtoV(U)
        return U, V
