#!/usr/bin/env/ python3

import numpy as np

def Output(V, nout):
    filename = './result/' + str(nout) + '.txt'
    np.savetxt(filename, V[0], delimiter=',')
    