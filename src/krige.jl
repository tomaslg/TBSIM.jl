#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function krige(datacoord,coord,model,cc,b,nugget,model_rotationmatrix,maxorder)
#---------------------------------------------
# Compute kriging weights at locations "coord"
#---------------------------------------------

# This program uses the following subroutine:
#     cova.m : compute covariance values for reduced distances

#------------------------------------------------------------------------------------------------------

# Definition of parameters
#-------------------------

    n = size(datacoord,1);  # number of data
    p = size(coord,1);      # number of points to estimate
    nst = size(model,1);    # number of nested structures

    j = Int(floor(gamma(maxorder+4)/gamma(maxorder+1+eps())/6 + 0.1)); # number of drift monomials
    if j>n
      println(" ");
      println("ERROR - Not enough data found in moving neighborhood. Conditioning kriging system is singular");
      println("        You should increase the search radius or the number of data per octant.");
      println(" ");
      return;
    end


# Calculation of the left covariance matrix K and right covariance matrix K0
#---------------------------------------------------------------------------

    x = vcat(datacoord,coord); # coordinates of data + locations to estimate

    k = zeros( convert(Dims,(n,n+p)) );

    for i = 1:nst

       # Calculation of matrix of reduced rotated distances
       R = model_rotationmatrix[:,:,i];
       t = x*R;
       t = t*t';
       # println("sqrt.(max.(0,-2*",t,"+","LinearAlgebra.diag(",t,")*ones(1,n+p)+ones(n+p,1)*LinearAlgebra.diag(",t,")'))")
       h = sqrt.(max.(0,-2*t+LinearAlgebra.diag(t)*ones(1,n+p)+ones(n+p,1)*LinearAlgebra.diag(t)'));
       h = h[1:n,:];

       # Evaluation of the current basic structure
       C = cova(model[i,1],h,b[i]);
       k = k+cc[i]*C;

    end

    k0 = k[:,(n+1):(n+p)]; # right member of the kriging system
    k = k[:,1:n];      # left member of the kriging system
    k = k+nugget*Matrix{Float64}(I, n, n);

    fl = [];
    fl0 = [];
    for i1 = 0:maxorder
      for i2 = 0:maxorder-i1
        for i3 = 0:maxorder-i1-i2
          fl = push!(fl,[(datacoord[jj,1]^i1)*(datacoord[jj,2]^i2)*(datacoord[jj,3]^i3) for jj=1:length(datacoord[:,3])]);
          fl0 = push!(fl0,[(coord[jj,1]^i1)*(coord[jj,2]^i2)*(coord[jj,3]^i3) for jj=1:length(coord[:,3])]);
        end
      end
    end
    fl=hcat(fl...);
    fl0=hcat(fl0...);

    if length(fl)>0
      k1=hcat(k,fl);
      if j>0
        k2=hcat(fl',zeros( convert(Dims, (j,j) ) ) );
      else
        k2=fl';
      end
    else
      k1=k;
      if j>0
        k2=zeros(convert(Dims,(j,j)));
      else
        k2=[];
      end
    end
    if length(k2)>0
      k = vcat(k1,k2);
    else
      k=k1;
    end
    if length(fl0)>0
      k0 = vcat(k0,fl0');
    end


    # Kriging weights
    #----------------
    # println("k=",k)
    # println("k0=",k0)
    # if det(k)==0
    #   weights = zeros(size(k0));
    # else
    weights = (k\k0);#inv(k)*k0
    # end
    # try
    # catch ee
    #   println("k0=",k0)
    #   println("k=",k)
    #   error(ee)
    # end
    weights = weights[1:n,1:p];
    return weights
end
