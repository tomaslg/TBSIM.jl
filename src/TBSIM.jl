module TBSIM

  export tableZY
  include("tbs.jl")
  # function tbsim_inc(simucoord=Float64[])
  #   return tbsim(simucoord)
  # end
  # export tbsim_inc
  export tbsim_spherical
  export tbsim_exponential
  export tbsim_gamma
  export tbsim_stable
  export tbsim_cubic
  export tbsim_Gaussian
  export tbsim_cardinal_sine
  export tbsim_J_Bessel
  export tbsim_K_Bessel
  export tbsim_generalized_Cauchy
  export tbsim_exponential_sine
  export tbsim_linear
  export tbsim_power
  function tbsim_spherical(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [1 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_spherical(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [1 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_exponential(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [2 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_exponential(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [2 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_gamma(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [3 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_gamma(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [3 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_stable(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [4 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_stable(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [4 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_cubic(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [5 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_cubic(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [5 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_Gaussian(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [6 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_Gaussian(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [6 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_cardinal_sine(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [7 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_cardinal_sine(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [7 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_J_Bessel(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [8 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_J_Bessel(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [8 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_K_Bessel(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [9 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_K_Bessel(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [9 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_generalized_Cauchy(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [10 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_generalized_Cauchy(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [10 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_exponential_sine(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [11 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_exponential_sine(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [11 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_linear(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [12 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_linear(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [12 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
  function tbsim_power(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
    return tbsim(
        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [13 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok
        )
  end
  function tbsim_power(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
    return tbsim(
        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
        x0,y0,z0,dx,dy,dz,nd,[-90,90],
        0.0,10.0,[1.0,5.0],1,0.0,[],
        [13 1.0 100 100 150 0 0 0 1 1000],
        9784498,[100 100 150] ,[0 0 0],
        1,4,3,ntok,true
        )
  end
#  function tbsim_mixed_power(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
#    return tbsim(
#        simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
#        1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
#        0.0,10.0,[1.0,5.0],1,0.0,[],
#        [14 1.0 100 100 150 0 0 0 1 1000],
#        9784498,[100 100 150] ,[0 0 0],
#        1,4,3,ntok
#        )
#  end
#  function tbsim_mixed_power(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
#    return tbsim(
#        Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
#        x0,y0,z0,dx,dy,dz,nd,[-90,90],
#        0.0,10.0,[1.0,5.0],1,0.0,[],
#        [14 1.0 100 100 150 0 0 0 1 1000],
#        9784498,[100 100 150] ,[0 0 0],
#        1,4,3,ntok,true
#        )
#  end
#  export tbsim_mixed_power
  # function tbsim_spline(nrealiz,simucoord,datacoord=[],ydata=[],ntok=5000)
  #   return tbsim(
  #       simucoord,datacoord,ydata,nrealiz,200,300,1,##nargin=2,
  #       1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
  #       0.0,10.0,[1.0,5.0],1,0.0,[],
  #       [15 0.45 100 100 150 0 0 0 2 1000],
  #       9784498,[100 100 150] ,[0 0 0],
  #       1,4,3,ntok
  #       )
  # end
  # function tbsim_spline(nrealiz,nx,ny,nz,x0,y0,z0,dx,dy,dz,datacoord=[],ydata=[],nd=Int16[1,1,1],ntok=5000)
  #   return tbsim(
  #       Float64[],[],[],nrealiz,nx,ny,nz,##nargin=2,
  #       x0,y0,z0,dx,dy,dz,nd,[-90,90],
  #       0.0,10.0,[1.0,5.0],1,0.0,[],
  #       [15 0.45 100 100 150 0 0 0 2 1000],
  #       9784498,[100 100 150] ,[0 0 0],
  #       1,4,3,ntok
  #       )
  # end
  # export tbsim_spline

  # α=2.5,simucoord=Float64[],datacoord=[],ydata=[],nrealiz = 10,return_a_mapping_to_grades=true,nargin=2,200,300,1,
  #     1.0,1.0,100.0,2.0,2.0,10.0,Int16[1,1,1],[-90,90],
  #     0.0,10.0,[1.0,5.0],1,0.0,[],
  #     [1 0.45 100 100 150 0 0 0 1 1000; 2 0.45 100 100 1000000000 0 0 0 1 1000],
  #     9784498,[100 100 150] ,[0 0 0],
  #     1,4,3,ntok = 5000
  #

end


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
