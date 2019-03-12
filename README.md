## TBSIM.jl

*TBSIM.jl* is a Julia package that implements a wide range of Gaussian field sampling methods.

USE: tbsim(simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail,model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok)

INPUT:
simucoord: coordinates of the locations to simulate (void if a grid simulation is required)
x0,y0,z0: if grid: minimum grid coordinates along x, y and z directions
nx,ny,nz:          number of grid nodes along x, y and z directions
dx,dy,dz:          grid meshes along x, y and z directions
nd: block discretization along x, y and z directions (1 \times 3 vector) ([1 1 1] for point-support simulation). The block size is given by the grid meshes dx,dy,dz.
datacoord: data coordinates (n * 3 matrix; void for non-conditional simulations)
ydata: Gaussian conditioning data (n * 1 vector; void for non-conditional simulations)
limits: trimming limits (inf and sup) for the Gaussian data (1 \times 2 vector)
tableZY: conversion table between original and Gaussian values (void if no transformation is required). The first column of the table contains the original values, the second column their normal scores
zmin,zmax: minimum and maximum values for the original variable
tail: parameters lambda and lambda prime for lower-tail and upper-tail extrapolations
model: covariance model for the Gaussian random field (nst \times 7 matrix, where nst is the number of nested structures). Each row refers to a nested structure and is codified as: [type, scale factors, angles]. There are three scale factors (along the rotated y, x and z axes) and three angles to define the coordinate rotation (azimuth, dip and plunge), see Deutsch and Journel, 1992, p. 25. Available types:
                    1: spherical
                    2: exponential
                    3: gamma (parameter b > 0)
                    4: stable (parameter b between 0 and 2)
                    5: cubic
                    6: Gaussian
                    7: cardinal sine
                    8: Bessel J (parameter b > 0)
                    9: Bessel K (parameter b > 0)
                   10: generalized Cauchy (parameter b)
                   11: exponential sine
                   12: linear generalized covariance
                   13: power generalized covariance (parameter b > 0)
                   14: mixed power generalized covariance (parameter b between 0 and 2)
                   15: spline generalized covariance (parameter b = even integer)
cc: sills or slopes of the nested structures (nst \times 1 vector)
b: third parameters for covariance definition (nst \times 1 vector) (used for covariance types 3,4,8,9,10 and 13)
nugget: nugget effect variance
nlines: number of lines to use for each nested structure (nst \times 1 vector)
nrealiz: number of realizations to draw
seed: seed number for generating random values
radius: maximum search radii along y, x and z (rotated system) for conditioning data. If radius = inf (infinite), a unique neighborhood is assumed and dual kriging is used
angles: angles for anisotropic search, according to the GSLIB conventions (Deutsch and Journel, 1992, p. 25)
octant: divide the neighborhood in octants? 1=yes, 0=no
ndata: number of conditioning data per octant (if octant=1) or in total (if octant=0)
name: name of output file, empty string will not write any.
nbdecimal: number of decimals for the values in output
header: create a GSLIB header in the output file? 1=yes, 0=no
ntok: maximum number of points to keep in memory and simulate simultaneously (optional). This number defines how many locations are projected onto the lines at each step of the simulation.


## Quick start

Some examples:

```julia
using TBSIM

#Sample Gaussian field using two nested structures for the correlation matrix, exponential and spherical.

# 

```
