#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def langevin(z, v, tau, eta, dt):
    # Implements numerical solution of the SDE
    # dz = v * dt
    # dv = -1/tau * v * dt + eta dW
    # by means of the Euler-Maruyama method
    # Equivalent to Eqs. (1.9a) and (1.9b) in Gillespie (1996),
    # Phys Rev E vol. 54 no. 2

    dW = np.random.normal(loc = 0, scale = np.sqrt(dt), size = z.shape)
    z_ = z + v*dt
    v_ = v - (1/tau)*v*dt + eta*dW
    return z_, v_


# Numerical parameters and system coefficients
Np   = 10000
Tmax = 1000
tau  = 10
eta  = 0.01
dt   = 0.1
Nt   = int(Tmax/dt) + 1
# Output interval
Nskip = 10

# Initial values
Z   = np.zeros(Np)
V   = np.zeros(Np)
var = np.zeros(Nt)

# File for output
txt = open('../output/variance_python.txt', 'w')

for i in range(1, Nt):
    Z, V = langevin(Z, V, tau, eta, dt)
    # write time and variance every Nskip timesteps
    if i % Nskip == 0:
        txt.write(f'{i*dt} {np.var(Z)}\n')

# Close output file
txt.close()
