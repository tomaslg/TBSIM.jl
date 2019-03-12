using TBSIM,Distributions, PyCall
#ENV["PYTHON"]="/home/tomas/anaconda3/bin/python3.6"
#"/home/tomas/anaconda3/bin/python3"
@pyimport matplotlib.pyplot as plt
matplotlibbackendpdf = pyimport("matplotlib.backends.backend_pdf")
@pyimport numpy as np
pp = matplotlibbackendpdf[:PdfPages]("simulation_plots.pdf")

uncond=[10 + 2,10 + 4,10 + 6 , 10 + 8, 10*3 + 2,10*3 + 4,10*3 + 6,10*3 + 8, 10*5 + 2,10*5 + 4,10*5 + 6,10*5 + 8
,10*7 + 2 ,10*7 + 4,10*7 + 6,10*7 + 8];
cond=setdiff(1:100,uncond);
##x=[ 1. 1. 1. ; 1. 2. 1. ; 1. 3. 1. ; 2. 1. 1. ; 2. 2. 1. ; 2. 3. 1. ; 3. 1. 1. ; 3. 2. 1. ; 3. 3. 1.];
x=Array{Float64}(undef,10^2,3);
for i=1:10
    for j=1:10
        x[(i-1)*10 + j,1]= i*1.;
        x[(i-1)*10 + j,2]= j*1.;
        x[(i-1)*10 + j,3]= 1.;
    end
end

## ;  1. 1. 2. ; 1. 2. 2. ; 1. 3. 2. ; 2. 1. 2. ; 2. 2. 2. ; 2. 3. 2. ; 3. 1. 2. ; 3. 2. 2. ; 3. 3. 2.];
y1=tbsim_spherical(1,x);
y2=tbsim_exponential(1,x);
y3=tbsim_gamma(1,x);
y4=tbsim_stable(1,x);
y5=tbsim_cubic(1,x);
y6=tbsim_Gaussian(1,x);
y7=tbsim_cardinal_sine(1,x);
y8=tbsim_J_Bessel(1,x);
y9=tbsim_K_Bessel(1,x);
y10=tbsim_generalized_Cauchy(1,x);
y11=tbsim_exponential_sine(1,x);
y12=tbsim_linear(1,x);
y13=tbsim_power(1,x);
#y14=tbsim_mixed_power(1,x[uncond,:]);

noescenarios=100
z1=mean(tbsim_spherical(noescenarios,x[cond,:],x[uncond,:],y1[uncond]),dims=2);
z2=mean(tbsim_exponential(noescenarios,x[cond,:],x[uncond,:],y2[uncond]),dims=2);
z3=mean(tbsim_gamma(noescenarios,x[cond,:],x[uncond,:],y3[uncond]),dims=2);
z4=mean(tbsim_stable(noescenarios,x[cond,:],x[uncond,:],y4[uncond]),dims=2);
z5=mean(tbsim_cubic(noescenarios,x[cond,:],x[uncond,:],y5[uncond]),dims=2);
z6=mean(tbsim_Gaussian(noescenarios,x[cond,:],x[uncond,:],y6[uncond]),dims=2);
z7=mean(tbsim_cardinal_sine(noescenarios,x[cond,:],x[uncond,:],y7[uncond]),dims=2);
z8=mean(tbsim_J_Bessel(noescenarios,x[cond,:],x[uncond,:],y8[uncond]),dims=2);
z9=mean(tbsim_K_Bessel(noescenarios,x[cond,:],x[uncond,:],y9[uncond]),dims=2);
z10=mean(tbsim_generalized_Cauchy(noescenarios,x[cond,:],x[uncond,:],y10[uncond]),dims=2);
z11=mean(tbsim_exponential_sine(noescenarios,x[cond,:],x[uncond,:],y11[uncond]),dims=2);
z12=mean(tbsim_linear(noescenarios,x[cond,:],x[uncond,:],y12[uncond]),dims=2);
z13=mean(tbsim_power(noescenarios,x[cond,:],x[uncond,:],y13[uncond]),dims=2);
#z14=mean(tbsim_mixed_power(10,x[cond,:],x[uncond,:],y14),dims=2);

Labels=["spherical","exponential","gamma","stable","cubic","Gaussian","cardinal_sine","J_Bessel","K_Bessel","generalized_Cauchy","exponential_sine","linear","power"]#,"mixed_power"]
Z=[z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13]#,z14]
Y=[y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13]#,y14]



for ℓ=1:length(Z)
  plt.clf()
  z=Y[ℓ]
  # fig, ax = plt.subplots()
  fig,(ax, ax2) = plt.subplots(2, figsize=(8, 16))
  f=unique(x[:,1])
  v=unique(x[:,2])
  M=Array{Float64}(undef,length(f),length(v))
  for i=1:length(f)
    for j=1:length(v)
      for k=1:length(z)
        if (x[:,1])[k]==f[i] && (x[:,2])[k]==v[j]
          M[i,j]=round( exp(z[k]/2.5),digits=3)
          break
        end
      end
    end
  end

  im = ax[:imshow](M)
  # We want to show all ticks...
  ax[:set_yticks](np.arange(length(f)))
  ax[:set_xticks](np.arange(length(v)))
  # ... and label them with the respective list entries
  ax[:set_yticklabels](f)
  ax[:set_xticklabels](v)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax[:get_xticklabels](), rotation=45, ha="right",
           rotation_mode="anchor")

  # Loop over data dimensions and create text annotations.
  for i=1:length(f)
      for j =1:length(v)
          text = ax[:text](j-1, i-1, M[i, j],
                         ha="center", va="center", color="w")
     end
   end

  ax[:set_title]("Number "*string(ℓ)*", "*Labels[ℓ]*", non-conditional")


  z=Z[ℓ]

  f=unique(x[:,1])
  v=unique(x[:,2])
  M=Array{Float64}(undef,length(f),length(v))

  for i=1:length(f)
    for j=1:length(v)
      for k=1:length(z)
        if (x[cond,1])[k]==f[i] && (x[cond,2])[k]==v[j]
          M[i,j]=round( z[k],digits=3)
          break
        end
        if k==length(z)
            for u in uncond
              if x[u,1]==f[i] && x[u,2]==v[j]
                M[i,j]=round( exp((Y[ℓ])[u]/2.5),digits=3)
                break
              end
            end
            # M[i,j]=0.
        end
      end
    end
  end

  im = ax2[:imshow](M)
  # We want to show all ticks...
  ax2[:set_yticks](np.arange(length(f)))
  ax2[:set_xticks](np.arange(length(v)))
  # ... and label them with the respective list entries
  ax2[:set_yticklabels](f)
  ax2[:set_xticklabels](v)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax2[:get_xticklabels](), rotation=45, ha="right",
           rotation_mode="anchor")

  # Loop over data dimensions and create text annotations.
  for i=1:length(f)
      for j =1:length(v)
          text = ax2[:text](j-1, i-1, M[i, j],
                         ha="center", va="center", color="w")
     end
   end

  ax2[:set_title]("Number "*string(ℓ)*", "*Labels[ℓ]*", conditional")


  fig[:tight_layout]()
  plt.savefig(pp, format="pdf")
end
pp[:close]()
