
# function interp1(xpt, ypt, x; method="linear", extrapvalue=nothing)
#
#     if extrapvalue == nothing
#         y = fill(0, size(x))#zeros(x)
#         idx = trues(size(x))
#         # println("idx=",idx)
#     else
#         y = extrapvalue*ones(size(x));
#         idx = (x .>= xpt[1]) .& (x .<= xpt[end])
#     end
#     if method == "linear"
#         intf = interpolate((xpt,), ypt, Gridded(Linear()))
#         y[idx] = Float64[intf(z,1) for z in x[idx]];
#
#     elseif method == "cubic"
#         itp = interpolate(ypt, BSpline(Cubic(Natural())), OnGrid())
#         intf = Interpolations.scale(itp, xpt)
#         y[idx] = [intf[xi] for xi in x[idx]]
#     end
#
#     return y
# end
#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function backtr(y,table,zmin,zmax,tail)

#----------------------------------------------------
# Back-transformation from Gaussian to original scale
#----------------------------------------------------

    p = size(table,1);
    m,n = size(y);
    z = zeros(m,n);


# Do not perform back-transformation if the conversion table is empty
#--------------------------------------------------------------------

    if isempty(table)
      z = y;
      return z;
    end

# Values in the lower tail (exponential extrapolation)
#-----------------------------------------------------

    z1 = table[1,1];
    y1 = table[1,2];
    I1 =  (LinearIndices(y))[findall( x->(x < y1), y)];


    if !isempty(I1)
      b0 = (z1-zmin)*exp.(-tail[1]*y1);
      z[I1] = zmin .+ b0*exp.(tail[1]*y[I1]);
    end

# Values in the upper tail (exponential extrapolation)
#-----------------------------------------------------

    zp = table[p,1];
    yp = table[p,2];
    # println("y=",y)

    I2 = (LinearIndices(y))[findall( x->(x > yp), y)];# find(y>yp);

    if !isempty(I2)
      bp = (zp-zmax)*exp(tail[2]*yp);
      z[I2] = zmax .+ bp*exp.(-tail[2]*y[I2]);
    end

# Within-class values
#--------------------

    I = sort(union(I1,I2));
    # try
    # catch ee
    #   if isa(ee,DimensionMismatch)
    #     println("I1,I2=",I1,I2);
    #   end
    #   error(ee);
    # end
    I3 = [ii for ii=1:(m*n)];
    deleteat!(I3,I);

    if !isempty(I3)
      # println("typeof(table)=",typeof(table));
      # open("t.txt", "w") do ff
      #   writedlm(ff,table);
      # end
      # table=readdlm("t.txt");
      intf = interpolate( (table[:,2],) , table[:,1] , Gridded(Linear()));
      z[I3] = Float64[intf(zz) for zz in y[I3] ]; # Table lookup
    end

    return z
end
