#!/usr/bin/env/ python3

import sys
import toml

from boundary import LeftFreeBoundary, RightFreeBoundary
from flux import HLL, HLLD
from grid import Grid
from initialcondition import InitialCondition
from output import Output
from reconstruction import Minmod, Ppm
from setdt import SetDt
from tvdrk import TVDRK2

if __name__ == '__main__':
    # set input parameters
    input_params = toml.load(open('input.toml'))
    # set order from reconstruction method
    if input_params['scheme']['rec'] == 'minmod' or input_params['scheme']['rec'] == 'mc':
        order = 2
    elif input_params['scheme']['rec'] == 'ppm' or input_params['scheme']['rec'] == 'ceno':
        order = 3
    else:
        sys.exit() 
    # make instance of grid
    grid = Grid()
    # set static grid
    cartesian = grid.static_cartesian(input_params['grid']['ix'], 
                                        order, 
                                        input_params['region']['x1'], 
                                        input_params['region']['x0'])
    # make instance of flux
    flux = HLLD()
    # make instance of reconstruct
    rec = Minmod(input_params['grid']['ix'])
    # make instance of boundary conditions
    xlbc = LeftFreeBoundary(order)
    xrbc = RightFreeBoundary(input_params['grid']['ix'], order)
    # make instance of TVDRK
    tvdrk = TVDRK2(flux, rec, xlbc, xrbc, 
                    input_params['grid']['ix'], 
                    cartesian['dx'], 
                    order)
    # make incetence of initial conditions
    ini = InitialCondition(input_params['grid']['ix'], order)
    # set initial condition
    U, V = ini.RJ2a()
    # set time step variables
    nstep = 0
    nout = 0
    t = 0.0
    tout = 0.0
    # main loop ----------
    while nstep < input_params['time']['nstop']:
        if t >= tout:
            print('time = '+str(t))
            Output(V, nout)
            tout += input_params['time']['tint']
            nout += 1
        if t >= input_params['time']['tstop']:
            print('time is over '+str(input_params['time']['tstop']))
            break
        nstep += 1
        print('nstep = '+str(nstep))
        dt = SetDt(V, cartesian['minlength'], input_params['time']['cfl'])
        print('dt = '+str(dt))
        U, V = tvdrk.time_step(U, V, dt)
        t += dt
    # ---------- main loop