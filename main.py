#!/usr/bin/env/ python3

import toml
from grid import Grid
from initialcondition import DW1aShockTube
from output import Output
from setdt import SetDt
from tvdrk import TVDRK

if __name__ == '__main__':
    # set input parameters
    input_params = toml.load(open('input.toml'))
    # make grid
    grid = Grid(input_params['grid']['ix'], 
                input_params['order']['order'], 
                input_params['region']['x1'], 
                input_params['region']['x0'])
    # set initial condition
    U, V = DW1aShockTube(grid['ixmax'], 
                            input_params['order']['order'])
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
        dt = SetDt(V, grid['minlength'], input_params['time']['cfl'])
        print('dt = '+str(dt))
        V = TVDRK(U, V, grid['ixmax'], grid['dx'], dt, input_params['order']['order'])
        t += dt
    # ---------- main loop