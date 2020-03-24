#!/usr/bin/env/ python3

import numpy as np

def UVtoF(U, V, b2, vb):
	F = np.zeros(U.shape)
	F[0] = U[1]
	F[1] = U[1] * V[1] + (V[7]+0.5*b2) - V[4] * V[4]
	F[2] = U[2] * V[1] - V[5] * V[4]
	F[3] = U[3] * V[1] - V[6] * V[4]
	F[4] = 0.0
	F[5] = V[5] * V[1] - V[2] * V[4]
	F[6] = V[6] * V[1] - V[3] * V[4]
	F[7] = (U[7]+V[7]+0.5*b2) * V[1] - vb * V[4]
	return F
