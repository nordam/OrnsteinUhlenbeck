# OrnsteinUhlenbeck
Just some code to look at some things about the Ornstein-Uhlenbeck process

## Fortran code
To compile the fortran code in the `src` folder, use for example

`gfortran -O3 -o ou ou.f90`

Running the resulting executable `ou` will write a time series of variance into the file `output/variance_fortran.txt`.

## Python code
To run the python code, simply run the script in the `src` folder

`python ou.py`

Running this script will write a time series of variance into the file `output/variance_python.txt`.

## Plotting
After running both the fortran and the python code, run the plotting script in the src folder

`python plot_both.py`

Running this script will read the two output files, and produce a plot of variance as a function of time.
