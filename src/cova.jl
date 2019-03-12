#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function cova(it,h,b)

#------------------------------------------
# Compute covariance for reduced distance h
#------------------------------------------

  epsilon = 1e-12;

  if (it < 2) # Spherical model
    C = 1 .- 1.5*min.(h,1) + 0.5*(min.(h,1).^3);
  elseif (it < 3) # Exponential model
    C = exp.(-h);
  elseif (it < 4) # Gamma model
    C = 1 ./ (1 .+ h).^b;
  elseif (it < 5) # Stable model
    C = exp.(-h.^b);
  elseif (it < 6) # Cubic model
    C = 1 .- 7*(min.(h,1).^2) + 35/4*(min.(h,1).^3) - 7/2*(min.(h,1).^5) + 3/4*(min.(h,1).^7);
  elseif (it < 7) # Gaussian model
    C = exp.(-h.^2);
  elseif (it < 8) # Cardinal sine model
    C = sin.(h)./(h .+ epsilon);
  elseif (it < 9) # J-Bessel model
    C = (h/2).^(-b)*gamma(b .+ 1).*besselj.(b,h .+ epsilon);
  elseif (it < 10) # K-Bessel model
    C = ((h/2).^b)*2/gamma(b).*besselk.(b,h .+ epsilon);
  elseif (it < 11) # Generalized Cauchy model
    C = 1 ./ (1 .+ h.^2).^b;
  elseif (it < 12) # Exponential sine model
    C = sin.(pi/2*exp.(-h));
  elseif (it < 13) # Linear model
    C = -h;
  elseif (it < 14) # Power model
    k = ceil(b/2);
    C = (-1).^k * h.^b;
  elseif (it < 15) # Mixed power model
    h = (h==1).*(h .+ epsilon) + (h!=1).*h;
    C = (1 .- h.^b)./log.(h);
  elseif (it < 16) # Spline model
    k = ceil(b/2);
    C = (-1).^(k+1) * h.^b .* log.(h .+ epsilon);
  else
    error("Unavailable covariance model");
  end

  return C
end
