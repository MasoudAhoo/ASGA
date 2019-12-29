ASGA package
=============


This matlab package provides a generic solver for Accelerated SubGradient Algorithm (ASGA) for solving structured nonsmooth convex problems of the form

  min_x f(x)+g(x), x in C

 for proper, convex, and lower semicontinuous f and g, and C is a simple constraint.

We provide a solver for solving the above-mentioned problem:
- ASGA1: single-subproble Accelerated SubGradient Algorithm;
- ASGA2: parameter-free single-subproble Accelerated SubGradient Algorithm;
- ASGA3: double-subprobleAccelerated SubGradient Algorithm;
- ASGA4: parameter-free double-subprobleAccelerated SubGradient Algorithm.

# Using the package

## Installation

All the functions in the package are given in the folder /MatlabCodes.

In order to use the function we recommend to execute the following command

```Matlab
addpath(genpath('.'))
```

if you are not working in the root folder of the package or replacing '.' by the location of the folder on your machine.


## Examples:

We recommend to look at the following files to see how to use the package:
* Demos/Demo_sparse_optimization.m: sparse recovery.

# References

[1] M. Ahookhosh, Accelerated first-order methods for large-scale convexoptimization: nearly optimal complexity under strongconvexity, Mathematical Methods of Operations Research, 89, 319--353 (2019).



