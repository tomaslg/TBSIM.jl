#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function setdual(model,coord,cc,b,datacoord,model_rotationmatrix,maxorder)

#-----------------------------------------------------------
# Calculate the right-hand side member of the kriging system
#-----------------------------------------------------------

# This program uses the following subroutine:
#     cova.m : compute covariance values for reduced distances


  if isempty(datacoord)
    return 0;
  end

  nst = size(model,1); # number of basic structures
  m0 = size(datacoord,1); # number of data
  k0 = zeros(1,m0);

  for i = 1:nst

    # Calculation of matrix of reduced rotated distances
    R = model_rotationmatrix[:,:,i];
    h = (ones(m0,1)*coord-datacoord)*R;
    h = h.^2;
    h = sqrt(sum(h'));

    # Evaluation of the current basic structure
    C = cova(model[i,1],h,b[i]);
    k0 = k0 + cc[i]*C;

  end

  fl0 = [];
  for i1 = 0:maxorder
    for i2 = 0:(maxorder-i1)
      for i3 = 0:(maxorder-i1-i2)
        fl0 = [fl0 (coord[1].^i1).*(coord[2].^i2).*(coord[3].^i3)];
      end
    end
  end
  k0 = vcat(k0,fl0');
  return k0;
end
