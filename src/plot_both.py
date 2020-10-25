#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt

times_py, variance_py   = np.loadtxt('../output/variance_python.txt').T
times_f90, variance_f90 = np.loadtxt('../output/variance_fortran.txt').T

fig = plt.figure(figsize = (7,4))

plt.plot(times_py,  variance_py,  label = 'Python result')
plt.plot(times_f90, variance_f90, '--', label = 'Fortran result')

plt.ylabel('Variance')
plt.xlabel('Time')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('../output/variances.png')
