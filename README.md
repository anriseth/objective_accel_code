# Objective acceleration and nonlinear GMRES for unconstrained optimization

This repository contains implementations of two optimization acceleration algorithms
discussed in *A. N. Riseth. Objective acceleration for unconstrained optimization, 2017*

The code contains modifications of the original Matlab code provided by Hans De Sterck on his
[website](http://math.uwaterloo.ca/amath-numerical-analysis-and-scientific-computing-group/hans-de-stercks-homepage/software-prof-hans-de-sterck/matlab-code-n-gmres-optimization).

**All references to N-GMRES-O in the code correspond to the O-ACCEL algorithm in _A. N. Riseth. Objective acceleration for unconstrained optimization_**

## Requirements
These files require the [Poblano toolbox](https://software.sandia.gov/trac/poblano) for MATLAB to run.

## List of relevant files

**`ngmres.m`** The N-GMRES algorithm as implemented by De Sterck

**`ngmres_o.m`** The Objective Acceleration algorithm by Riseth, following the same structure as `ngmres.o`.

**`ngmres_test_general.m`** Provides convergence plots for the different algorithms.
To show convergence plots for an instance of Problem A, with n=200, call
``` matlab
ngmres_test_general(0,1,200,400,true)
```

**`runme_writestats.m`** Runs 1000 instances of each test problem and writes the statistics to file in the directory `data`.