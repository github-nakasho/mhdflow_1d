#!/usr/bin/env/ python3

import toml
from grid import Grid
from initialcondition import DW1aShockTube

if __name__ == '__main__':
    # set input parameters
    input_params = toml.load(open('input.toml'))
    # make gtid
    grid = Grid(input_params['grid']['ix'], 
                input_params['order']['order'], 
                input_params['region']['x1'], 
                input_params['region']['x0'])
    # set initial condition
    U, V = DW1aShockTube(grid['ixmax'])
    '''
    # set time step variables
    nstep = 0
    t = 0.0
    nout = 0
    # main loop ----------
    while nstep < nstop + 1:
        nstep += 1
        print('nstep = '+str(nstep))
        dt = SetDt(U, V)
        print('dt = '+str(dt))
        TVDRK()
        t += dt
        if t > tout:
            print('time = '+str(t))
			tout += tint
			nout += 1
		}
		if t > tstop:
			print('time is over '+str(tstop))
    # ---------- main loop
    '''