#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------
function fix(t)
    if t>=0
        return floor(t)
    else
        return ceil(t)
    end
end
function vdc(nlines,nrealiz,seed)

#------------------------------------------------------------
# Generation of equidistributed lines over the unit 3D sphere
# according to a van der Corput sequence
#------------------------------------------------------------
#
# USE: lines = vdc(nlines,nrealiz,seed)
#
# INPUT:
#   nlines: number of lines to generate for each realization
#   nrealiz: number of realizations
#   seed: seed for generation of random numbers
#
# OUTPUT:
#   lines: n*3 matrix representing the lines

# This program uses the following subroutine:
#     setrot.m: set up matrix for rotation and reduction of coordinates
#HERE I AM
    lines = zeros(nlines*nrealiz,3);
    rng = Random.MersenneTwister(seed);
    # rand('state',seed);
    # seed2 = abs(rand(Int32));#ceil(1e7*rand(rng));
    # rng = Random.MersenneTwister(seed2);
    # randn(rng);
    # randn(rng, Float64)

    i = Float32[j for j=1:nlines];

    # binary decomposition of i
    j = i;
    u = 0;
    p = 0;
    while (maximum(j) > 0)
      p = p+1;
      t = fix.(j/2);
      u = u .+ 2*(j/2 - t) ./ (2 .^ p);
      j = t;
    end

    # ternary decomposition of i
    j = i;
    v = 0;
    p = 0;
    while (maximum(j)>0)
      p = p+1;
      t = fix.(j/3);
      v = v .+ 3*(j/3 - t)./(3 .^ p);
      j = t;
    end

    # directing vector of the i-th line
    x  = [cos.(2*pi*u).*sqrt.(1 .- v.*v) sin.(2*pi*u).*sqrt.(1 .- v.*v) v];

    # random rotation
    for k = 1:nrealiz
      angles = 360*rand(rng,1,3);
      R = setrot(hcat([1 1 1 1] ,angles),1);
      lines[(k-1)*nlines+1:k*nlines,:] = x*R;
    end
    return lines
end
