#!/usr/bin/env/ python3

import sys
import toml

from boundary import LeftFreeBoundary, RightFreeBoundary
from flux import HLL, HLLD
from grid import Grid
from initialcondition import InitialCondition
from output import Output
from reconstruction import Minmod, MC, VanLeer, Ppm, CENO, MP5
from setdt import SetDt
from tvdrk import TVDRK2, TVDRK3

if __name__ == '__main__':
    # set input parameters
    input_params = toml.load(open('input.toml'))
    # set order from reconstruction method
    if input_params['scheme']['rec'] in ['minmod', 'mc', 'vanleer']:
        order = 2
    elif input_params['scheme']['rec'] in ['ppm', 'ceno']:
        order = 3
    elif input_params['scheme']['rec'] in ['mp5']:
        order = 5
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
    if input_params['scheme']['flux'] == 'HLL':
        flux = HLL()
    elif input_params['scheme']['flux'] == 'HLLD':
        flux = HLLD()
    else:
        sys.exit()
    # make instance of reconstruct
    rec = MP5(input_params['grid']['ix'], order)
    # make instance of boundary conditions
    xlbc = LeftFreeBoundary(order)
    xrbc = RightFreeBoundary(input_params['grid']['ix'], order)
    # make instance of TVDRK
    if input_params['scheme']['tvdrk'] == 2:
        tvdrk = TVDRK2(flux, rec, xlbc, xrbc, 
                        input_params['grid']['ix'], 
                        cartesian['dx'], 
                        order)
    elif input_params['scheme']['tvdrk'] == 3:
        tvdrk = TVDRK3(flux, rec, xlbc, xrbc, 
                        input_params['grid']['ix'], 
                        cartesian['dx'], 
                        order)
    # make incetence of initial conditions
    ini = InitialCondition(cartesian['x'], input_params['grid']['ix'], order)
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
    