#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function dual(datacoord,ydata,model,cc,b,nugget,model_rotationmatrix,maxorder)

#-----------------------------
# Compute dual kriging weights
#-----------------------------

# This program uses the following subroutine:
#     cova.m : compute covariance values for reduced distances

#------------------------------------------------------------------------------------------------------

# Definition of parameters
#-------------------------
    n = size(datacoord,1);  # number of data
    nst = size(model,1);    # number of nested structures
    j = floor(gamma(maxorder+4)/gamma(maxorder+1+eps)/6 + 0.1); # number of drift monomials
    if j>n
      println(" ");
      println("ERROR - Not enough data found in unique neighborhood. Conditioning kriging system is singular");
      println("        You should lower the exponent of the power model.");
      println(" ");
      return ;
    end


    # Calculation of the left covariance matrix K and right covariance matrix K0
    #---------------------------------------------------------------------------

    k = zeros(n,n);

    for i = 1:nst

       # Calculation of matrix of reduced rotated distances
       R = model_rotationmatrix[:,:,i];
       t = datacoord*R;
       t = t*t';
       h = sqrt(max(0,-2*t+LinearAlgebra.diag(t)*ones(1,n)+ones(n,1)*LinearAlgebra.diag(t)'));

       # Evaluation of the current basic structure
       C = cova(model[i,1],h,b[i]);
       k = k+cc[i]*C;

    end

    k = k+nugget*LinearAlgebra.eye(n);

    fl = [];
    for i1 = 0:maxorder
      for i2 = 0:(maxorder-i1)
        for i3 = 0:(maxorder-i1-i2)
          fl = hcat(fl, (datacoord[:,1].^i1).*(datacoord[:,2].^i2).*(datacoord[:,3].^i3) );
        end
      end
    end
    k = vcat(hcat(k,fl),hcat(fl',zeros(j,j)));
    ydata = vcat(ydata,zeros(j,1));


    # Kriging weights
    #----------------

    weights = ydata\k;# =inv(k)*ydata
    return weights
end
