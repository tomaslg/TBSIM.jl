using TBSIM

uncond=[2,5,8];
cond=[1,3,4,6,7,9];
x=[ 1. 1. 1. ; 1. 2. 1. ; 1. 3. 1. ; 2. 1. 1. ; 2. 2. 1. ; 2. 3. 1. ; 3. 1. 1. ; 3. 2. 1. ; 3. 3. 1.];
## ;  1. 1. 2. ; 1. 2. 2. ; 1. 3. 2. ; 2. 1. 2. ; 2. 2. 2. ; 2. 3. 2. ; 3. 1. 2. ; 3. 2. 2. ; 3. 3. 2.];
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
#y14=tbsim_mixed_power(1,x[uncond,:]);


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
#z14=tbsim_mixed_power(1,x[cond,:],x[uncond,:],y14);


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
#w14,pos14=tbsim_mixed_power(1,20,20,3,0.,0.,0.,5.,5.,5.)

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
#v14,Pos14=tbsim_mixed_power(1,19,19,3,2.,2.,0.,5.,5.,5.,pos14[uncon,:],w14[uncon])

# # #@enter
# # outputformat,datacoord,extreme_coord,ydata,simucoord,COORDS,RR=tbsim()
# # Λ=readdlm("nscore.out");
# # datacoord=Λ[:,1:3];
# # ydata=Λ[:,4];
# wdir=dirname(dirname(pwd()))*"/data/";
# lws=readdlm(wdir*"scenarios/laws.scen");
# ydata=lws[1:12000];
# ξnew=lws[12001:end];
# pos=readdlm(wdir*"scenarios/Grilla.dat");
# datacoord=pos[1:12000,:];
# simucoord=pos[12001:end,:];
#
# Time=time();
# outputformat=tbsim(2.5,simucoord,datacoord,ydata);
# outputformat1=tbsim(2.5,simucoord);
# println("time[sec]=",time()-Time);
# # println("outputformat1=",mean(outputformat,1))
# # println("outputformat2=",mean(outputformat,2))
# # println("std1=",std(outputformat,1))
# # println("st22=",std(outputformat,2))
# # writedlm("outputformat.txt",outputformat,";");
# #
# # coord=vcat(COORDS...);
# #
# # ENV["PYTHON"]="python3.6"
# # using PyCall
# #
# # @pyimport matplotlib.pyplot as plt
# # matplotlibbackendpdf = pyimport("matplotlib.backends.backend_pdf")
# # pp = matplotlibbackendpdf[:PdfPages]("allrealizationshm.pdf")
# #
# # for ξ=1:10
# #   hm=Array{Float64}(300,200)
# #   k=0#300*200
# #   for i=1:300
# #     for j=1:200
# #       k=k+1;
# #       hm[i,j]=outputformat[k,ξ]#mean(outputformat[k,:])
# #     end
# #   end
# # 	plt.clf()
# #   plt.imshow(hm, cmap="hot", interpolation="nearest")
# #   plt.savefig(pp, format="pdf")
# # end
# # hm=Array{Float64}(300,200)
# # k=0#300*200
# # for i=1:300
# #   for j=1:200
# #     k=k+1;
# #     hm[i,j]=mean(outputformat[k,:])
# #   end
# # end
# # plt.clf()
# # plt.imshow(hm, cmap="hot", interpolation="nearest")
# # # plt.savefig(pp, format="pdf")
# # # k=0#300*200
# # # for i=1:300
# # #   for j=1:200
# # #     k=k+1;
# # #     hm[i,j]=mean(RR[k])
# # #   end
# # # end
# # # plt.clf()
# # # plt.imshow(hm, cmap="hot", interpolation="nearest")
# # plt.savefig(pp, format="pdf")
# # pp[:close]()
