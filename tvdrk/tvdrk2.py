#!/usr/bin/env/ python3

import numpy as np

from convert import Convert
from source import Source


class TVDRK2:
    def __init__(self, flux, rec, xlbc, xrbc, ix, dx, order):
        self.flux = flux
        self.rec = rec
        self.xlbc = xlbc
        self.xrbc = xrbc
        self.ix = ix
        self.dx = dx
        self.order = order
        self.ixmax = ix + 2 * (order-1)
        self.conv = Convert()
        
    def half_step(self, U, V, dt):
        ix = self.ix
        dx = self.dx
        order = self.order
        U1 = np.zeros(U.shape)
        Vl, Vr = self.rec.make_rec(V)
        Vl = self.xlbc.set_vl(Vl, Vr)
        Vr = self.xrbc.set_vr(Vl, Vr)
        F = self.flux.make_flux(Vl, Vr, ix, order)
        S = Source(U)
        for m in range(8):
            for i in range(order-1, ix+(order-1)):
                U1[m][i] = U[m][i] - dt / dx * (F[m][i]-F[m][i-1]) + dt * S[m][i]
        U1 = self.xlbc.set_u(U1)
        U1 = self.xrbc.set_u(U1)
        V = self.conv.UtoV(U1)
        return U1, V

    def full_step(self, U, U1, V, dt):
        ix = self.ix
        dx = self.dx
        order = self.order
        Vl, Vr = self.rec.make_rec(V)
        Vl = self.xlbc.set_vl(Vl, Vr)
        Vr = self.xrbc.set_vr(Vl, Vr)
        F = self.flux.make_flux(Vl, Vr, ix, order)
        S = Source(U1)
        for m in range(8):
            for i in range(order-1, ix+(order-1)):
                U[m][i] = 0.5 * (U[m][i]+U1[m][i]-dt/dx*(F[m][i]-F[m][i-1])+dt*S[m][i])
        U = self.xlbc.set_u(U)
        U = self.xrbc.set_u(U)
        V = self.conv.UtoV(U)
        return U, V

    def time_step(self, U, V, dt):
        U1, V = self.half_step(U, V, dt)
        U, V = self.full_step(U, U1, V, dt)
        return U, V
