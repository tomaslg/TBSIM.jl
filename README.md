[![Build Status](https://travis-ci.org/tomaslg/TBSIM.jl.svg?branch=master)](https://travis-ci.org/tomaslg/TBSIM.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/covxo242g63naqgm?svg=true)](https://ci.appveyor.com/project/tomaslg/tbsim-jl)
[![Coverage Status](https://coveralls.io/repos/github/tomaslg/TBSIM.jl/badge.svg?branch=master)](https://coveralls.io/github/tomaslg/TBSIM.jl?branch=master)

##TBSIM.jl

*TBSIM.jl* is a Julia package that implements a wide range of Gaussian field sampling methods.

[Paper](https://www.sciencedirect.com/science/article/pii/S0098300406000549)

USE: tbsim(simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail,model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok)

#INPUT:
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

uncond=[2,5,8];#index of points that'll affect the support of the conditional simulation
cond=[1,3,4,6,7,9];#index of points that will be sampled a-posteriori
x=[ 1. 1. 1. ; 1. 2. 1. ; 1. 3. 1. ; 2. 1. 1. ; 2. 2. 1. ; 2. 3. 1. ; 3. 1. 1. ; 3. 2. 1. ; 3. 3. 1.];
#Matrix of points
```
The following calls return an unconditional simulation at the 3D points in matrix x using different variogram structures.
```julia
y1=tbsim_spherical(1,x[uncond,:]);
y2=tbsim_exponential(1,x[uncond,:]);
y3=tbsim_gamma(1,x[uncond,:]);
y4=tbsim_stable(1,x[uncond,:]);
y5=tbsim_cubic(1,x[uncond,:]);
y6=tbsim_Gaussian(1,x[uncond,:]);
y7=tbsim_cardinal_sine(1,x[uncond,:]);
y8=tbsim_J_Bessel(1,x[uncond,:]);
y9=tbsim_K_Bessel(1,x[uncond,:]);
y10=tbsim_generalized_Cauchy(1,x[uncond,:]);
y11=tbsim_exponential_sine(1,x[uncond,:]);
y12=tbsim_linear(1,x[uncond,:]);
y13=tbsim_power(1,x[uncond,:]);
```
Then you can sample an scenario of a conditional simulations like this:
```julia
z1=tbsim_spherical(1,x[cond,:],x[uncond,:],y1);
z2=tbsim_exponential(1,x[cond,:],x[uncond,:],y2);
z3=tbsim_gamma(1,x[cond,:],x[uncond,:],y3);
z4=tbsim_stable(1,x[cond,:],x[uncond,:],y4);
z5=tbsim_cubic(1,x[cond,:],x[uncond,:],y5);
z6=tbsim_Gaussian(1,x[cond,:],x[uncond,:],y6);
z7=tbsim_cardinal_sine(1,x[cond,:],x[uncond,:],y7);
z8=tbsim_J_Bessel(1,x[cond,:],x[uncond,:],y8);
z9=tbsim_K_Bessel(1,x[cond,:],x[uncond,:],y9);
z10=tbsim_generalized_Cauchy(1,x[cond,:],x[uncond,:],y10);
z11=tbsim_exponential_sine(1,x[cond,:],x[uncond,:],y11);
z12=tbsim_linear(1,x[cond,:],x[uncond,:],y12);
z13=tbsim_power(1,x[cond,:],x[uncond,:],y13);


```

If you want to generate simulations at grid points you may call the functions like as follows:
```julia

w1,pos1=tbsim_spherical(1,20,20,3,0.,0.,0.,5.,5.,5.)
w2,pos2=tbsim_exponential(1,20,20,3,0.,0.,0.,5.,5.,5.)
w3,pos3=tbsim_gamma(1,20,20,3,0.,0.,0.,5.,5.,5.)
w4,pos4=tbsim_stable(1,20,20,3,0.,0.,0.,5.,5.,5.)
w5,pos5=tbsim_cubic(1,20,20,3,0.,0.,0.,5.,5.,5.)
w6,pos6=tbsim_Gaussian(1,20,20,3,0.,0.,0.,5.,5.,5.)
w7,pos7=tbsim_cardinal_sine(1,20,20,3,0.,0.,0.,5.,5.,5.)
w8,pos8=tbsim_J_Bessel(1,20,20,3,0.,0.,0.,5.,5.,5.)
w9,pos9=tbsim_K_Bessel(1,20,20,3,0.,0.,0.,5.,5.,5.)
w10,pos10=tbsim_generalized_Cauchy(1,20,20,3,0.,0.,0.,5.,5.,5.)
w11,pos11=tbsim_exponential_sine(1,20,20,3,0.,0.,0.,5.,5.,5.)
w12,pos12=tbsim_linear(1,20,20,3,0.,0.,0.,5.,5.,5.)
w13,pos13=tbsim_power(1,20,20,3,0.,0.,0.,5.,5.,5.)

```

The call of the function recieves (number of scenarios to sample(1), number of discretizations for the grid on the x(20),y(20),z(3) axis, origen(0.,0.,0.), step (dx,dy,dz) on each axis (5.,5.,5.)), and returns a Float64 with the values of the simulation and an Array{Float64,2} of 3D points. Similarly as in the previous calls, conditional simulations can be sampled like this:

```julia

uncon=[4*i for i=1:Int(floor(length(w13)/4))]

v1,Pos1=tbsim_spherical(1,19,19,3,2.,2.,0.,5.,5.,5.,pos1[uncon,:],w1[uncon])
v2,Pos2=tbsim_exponential(1,19,19,3,2.,2.,0.,5.,5.,5.,pos2[uncon,:],w2[uncon])
v3,Pos3=tbsim_gamma(1,19,19,3,2.,2.,0.,5.,5.,5.,pos3[uncon,:],w3[uncon])
v4,Pos4=tbsim_stable(1,19,19,3,2.,2.,0.,5.,5.,5.,pos4[uncon,:],w4[uncon])
v5,Pos5=tbsim_cubic(1,19,19,3,2.,2.,0.,5.,5.,5.,pos5[uncon,:],w5[uncon])
v6,Pos6=tbsim_Gaussian(1,19,19,3,2.,2.,0.,5.,5.,5.,pos6[uncon,:],w6[uncon])
v7,Pos7=tbsim_cardinal_sine(1,19,19,3,2.,2.,0.,5.,5.,5.,pos7[uncon,:],w7[uncon])
v8,Pos8=tbsim_J_Bessel(1,19,19,3,2.,2.,0.,5.,5.,5.,pos8[uncon,:],w8[uncon])
v9,Pos9=tbsim_K_Bessel(1,19,19,3,2.,2.,0.,5.,5.,5.,pos9[uncon,:],w9[uncon])
v10,Pos10=tbsim_generalized_Cauchy(1,19,19,3,2.,2.,0.,5.,5.,5.,pos10[uncon,:],w10[uncon])
v11,Pos11=tbsim_exponential_sine(1,19,19,3,2.,2.,0.,5.,5.,5.,pos11[uncon,:],w11[uncon])
v12,Pos12=tbsim_linear(1,19,19,3,2.,2.,0.,5.,5.,5.,pos12[uncon,:],w12[uncon])
v13,Pos13=tbsim_power(1,19,19,3,2.,2.,0.,5.,5.,5.,pos13[uncon,:],w13[uncon])

```
