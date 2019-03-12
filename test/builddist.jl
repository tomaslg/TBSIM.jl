#Calculate first μ,Σ from data

#include("initparamsDMP.jl")

@everywhere cff=10
@everywhere crf=0.25
@everywhere cm=2.5
@everywhere recovery=0.85
@everywhere copperprice=2.1
@everywhere tonpound=2204.62

function buildlawscen(scenarioslaws,Ω,NoEscenarios,nblocks)
    scen=zeros(nblocks,NoEscenarios)
    for i =1:nblocks
        for ξ in 1:NoEscenarios
            scen[i,ξ] = ((scenarioslaws[i,2*ξ]/Ω[i]) + cff+cm)/(recovery*(copperprice-crf)*tonpound)
        end
    end
    return scen
end

@everywhere function f(ar)
    row1,row2,E1,E2=ar
    return dot(row1,row2)/(length(row1)-1) - length(row1)*E1*E2/(length(row1)-1)
end

function builddistcluster(scen,NoEscenarios,𝝹)
    μ = sum(scen[:,ξ] for ξ = 1:NoEscenarios)/NoEscenarios
    # args=[]
    a=[]
    for i =1:𝝹
        for j = i:𝝹
            push!(a,f([ [ scen[i,k] for k=1:NoEscenarios ] ,  [ scen[j,k] for k=1:NoEscenarios ] ,μ[i],μ[j]]))
        end
    end
    # a = SharedArray{Float32}(Int(𝝹*(𝝹+1)/2))
    # #
    # # a=zeros(Int(𝝹*(𝝹+1)/2))
    # @parallel for k = 1:Int(𝝹*(𝝹+1)/2)
    #     a[k] = f(args[k])
    # end

    # 𝐳=pmap(f,args)
    # a=fetch(𝐳)

    k=1
    R=zeros(Float32,𝝹,𝝹)
    for i =1:𝝹
        for j = i:𝝹
            R[i,j]=R[j,i]=(a[k])
            k=k+1
        end
    end
    return μ,R
end

@everywhere function linearTprc(lawscc,Ω,blockcluster,NoEscenarios,bci)
    scenarioslaws=zeros(length(Ω),2*NoEscenarios)
    for ξ =1:NoEscenarios
        for i in bci
            scenarioslaws[i,2*ξ - 1] =  round(-cm*Ω[i],digits=6)
            scenarioslaws[i,2*ξ] = round((max(lawscc[i,ξ],0.)*tonpound*(copperprice-crf)*recovery - cff -cm)*Ω[i],digits=6)
        end
    end
    return scenarioslaws
end


function builddist(scen,Ξ,nblocks,clusterblocks)
    L= zeros(length(clusterblocks),length(Ξ))
    Ωc=Float64[]
    sc= zeros(scen)#[scenarioslaws[:,2*ξ] for ξ in Ξ ]
    for (k,c) in enumerate(clusterblocks)
        for ξ in 1:length(Ξ)
            L[k,ξ]=sum(scen[i,ξ] for i in c)/length(c)
            for i in c
                sc[i,ξ] =scen[i,ξ] - L[k,ξ]
            end
        end
    end
    μc,Σc = builddistcluster(L,length(Ξ),length(clusterblocks))
    μb = zeros(nblocks)
    Σb = zeros(nblocks)
    for i=1:nblocks
        μb[i] = sum(sc[i,ξ] for ξ in Ξ)/length(Ξ)
        Σb[i] = sum(sc[i,ξ]^2 - length(Ξ) * μb[i]^2 for ξ in Ξ)/(length(Ξ)-1)
    end
    return μc,Σc,μb,Σb
end


#Update distribution of N(μ,Σ), given by:
#   #Σ1|2= Σ11 - Σ12 Σ22^{-1} Σ21
#   #μ1|2 = μ1 + Σ12 Σ22^{-1} (ξ2-μ2)


#
# function Σ12invR22ξ2μ2(μ,lawsccξ,C,EC,R)
#     sd=setdiff(C,EC)
#     ξ2μ2 = Float64[]
#     ξ2μ2 = [lawsccξ[i]-μ[i] for i in EC]
#     L= transpose(R[sd,EC])\R[EC,EC]
#     return  L, L * ξ2μ2
# end

# function updatedistribution(μ,C,EC,R,lawsccξ)
# #μ mean, C set of clusters, EC set of extracted clusters before τ-1, R variance-covariance Matrix,
#     Ƶ=copy(R)
#     Σ12iR22,Σ12iR22ξ2μ2=Σ12invR22ξ2μ2(μ,lawsccξ,C,EC,Ƶ)
#     sd=setdiff(C,EC)
#     μ1= Float64[]
#     j=1
#     for i in C
#         if i in sd
#             push!(μ1,μ[i]+Σ12iR22ξ2μ2[j])
#             j=j+1
#         else
#             push!(μ1,μ[i])
#         end
#     end
#     Σ = Ƶ[sd,sd] - Σ12iR22 * Ƶ[EC,sd]
#     ki=1
#     for i in sd
#         kj=1
#         for j in sd
#             Ƶ[i,j]=Ƶ[j,i]= Σ[ki,kj]
#             kj=kj+1
#         end
#         ki=ki+1
#     end
#     return μ1,Ƶ
# end

# Given a multivariate distribution p(x),
# returns new distribution p(x | xᵢ=aᵢ), where
# "a" is a vector of values to condition on, and
# "d2" is a vector of indices providing the
# dimensions (subscript i) in increasing order.
function updatedistribution1(
    joint::MultivariateDistribution,
    a::Vector{Float64},
    d2::Vector{Int32}
    )
    if length(d2)==0 return joint.μ,Matrix(joint.Σ) end
    if length(d2) != length(a)
        error("list of conditioned values and dimension indices bet the same length")
    end
    max_dim = length(joint.μ)
    if minimum(d2) < 1 || maximum(d2) > max_dim
        error("dimension indices must be between 1 and dimension of joint distribution")
    end

    # d1 = remaining dimensions
    # d2 = removed/conditioned dimensions
    d2 = sort(d2)
    d1 = setdiff(1:max_dim,d2)

    # covariance matrix blocks
    Σ = Matrix(joint.Σ);#PDMats.full(joint.Σ)
    Σ11 = Hermitian(Σ[d1,d1])
    Σ12 = Σ[d1,d2]
    Σ22inv = Hermitian(inv(Σ[d2,d2]))

    μ_new = joint.μ[d1] + (Σ12 * Σ22inv * (a - joint.μ[d2]))
    Σ_new = Σ11 - X_invA_Xt(PDMat(Σ[d2,d2]), Σ12)# (Σ12 * Σ22inv * transpose(Σ12))
    return μ_new,Σ_new
end

@everywhere function updatedistribution2(dhb,bb,μb,Σb,surface_extracted,surface_clusterblocks,ξ,DH,blockcluster,EC)
#
    μprime=copy(μb)
    Σprime=copy(Σb)
    surface_blocks=length(surface_extracted)==0 ? [] : [j for i in surface_extracted if i!=0 for j in surface_clusterblocks[i]]
    # for b in surface_blocks
    #     μprime[b]=ξ[b]
    #     Σprime[b]=0.
    # end
    for b =1:length(bb)
        if !(blockcluster[b] in EC)
            Γ=intersect(bb[b],surface_blocks)
            if length(Γ)>0
                μprime[b]= mean(vcat(DH[dhb[b]],ξ[Γ]))
                Σprime[b]=std(vcat(DH[dhb[b]],ξ[Γ]))
                # sqrt(( sum(DH[h]^2 for h in dhb[b]) + sum(ξ[h]^2 for h in Γ) -
                #             (sum(DH[h] for h in dhb[b]) + sum(ξ[h] for h in Γ))^2 / (length(dhb[b]) + length(Γ))   )/(length(dhb[b]) + length(Γ)-1))
            else
                μprime[b]=mean(DH[dhb[b]])
                Σprime[b]= length(dhb[b])<=1 ? 0.02 : std(DH[dhb[b]])
                # sqrt(( sum(DH[h]^2 for h in dhb[b])  -
                # (sum(DH[h] for h in dhb[b]))^2 / length(dhb[b]) )/(length(dhb[b]) - 1))
            end
        end
    end
    return μprime,Σprime
end

function choleskyfactor(R)#R real valued
    L = zeros(R)
    for i =1:length(R[:,1])
        L[i,1]= i!=1 ? (1/L[1,1])*R[i,1] : sqrt(R[1,1])
    end
    for j =2:(length(R[1,:])-1)
        L[j,j]= sqrt( R[j,j] - sum( L[j,k]^2 for k=1:(j-1)) )
        for i =(j+1):length(R[:,j])
            L[i,j]= (1/L[j,j])*(R[i,j] - sum( L[i,k]*L[j,k] for k=1:(j-1)) )
        end
    end
    L[end,end]=sqrt( R[end,end] - sum( L[end,k]^2 for k=1:(length(R[1,:])-1)) )
    return L
end

# function sample(μc,R,Ξ,C,EC,clusterblocks)
#     sd=setdiff(C,EC)
#     # MVN = MvNormal(μc[sd],R[sd,sd])
#     # D=rand(MVN,length(Ξ))
#     lawscc = zeros(length(μc),length(Ξ))
#     lawscb = zeros(length(μb),length(Ξ))
#     Ch = (cholfact( Symmetric(R[sd,sd]) ))[:U]#choleskyfactor(R[sd,sd])#
#     for ξ in Ξ
#         rd=randn(length(sd))#rd[1]=randn()
#         D= μc[sd] + Ch*rd
#         for (j,i) in enumerate(sd)
#             lawscc[i,ξ]= D[j]
#             for 𝐣 in clusterblocks[i]
#                 lawscb[𝐣,ξ] = μb[𝐣] #(μb[𝐣]+μc[j])*lawscc[j,ξ]/μc[j] - lawscc[j,ξ] #+ randn()*Σb[𝐣] μb[𝐣]#
#                 #*recovery*(copperprice-crf)*tonpound
#                 #- cff - cm)
#                 #*sum(Ω[𝐢] for 𝐢 in clusterblocks[i])/length(clusterblocks[i])
#                 #+(μb[𝐣] + randn()*Σb[𝐣] )*Ω[𝐣]
#             end
#         end
#     end
#     # printsolution(wx,bench,phase,[μc for τ=1:T])
#     # μ=[sum(lawscc[i,:]) for i = 1:length(lawscc[:,1])]/length(lawscc[1,:])
#     # printsolution(wx,bench,phase,[μ for τ=1:T])
#     return lawscc,lawscb
# end

function updatesamples1(lawsccξ,Ξ,C,EC#=lawsccξ,Ξ,C,EC=Float64[j in ECadp ? newξ[j] : lawsccadp[j,ξ] for j in C ],1:InSample,C,[j for j in union(ECadp,ECτ)]=#,
    MVClusters,meanscenarioinsampleadp)
    μ_new,Σ_new=updatedistribution1(MVClusters,lawsccξ[EC],EC)
    MVClustersnew = MvNormal(μ_new,Σ_new)
    # Σ12iR22=Σ12invR22(C,EC,ECτ,R)
    if meanscenarioinsampleadp
        lawscc=mean(MVClustersnew)#μnew#,lawscb,μb
    else
        lawscc=rand(MVClustersnew,length(Ξ)) #,lawscb=sample(μnew,Σnew,μb,Σb,Ξ,C,EC,clusterblocks)
    end
    sample=Array{Float32}(undef,length(C),meanscenarioinsampleadp ? 1 : length(Ξ))
    for ξ in (meanscenarioinsampleadp ? 1 : Ξ)
        sample[:,ξ]=mean(MVClusters)
    end
    j=0
    for i in C
        if !(i in EC)
            j=j+1
            sample[i,:]=lawscc[j,:]
        end
    end
    # for ii in ECτ
    #     lawscc[ii]=lawsccξ[ii]#μc[ii]
    # end
    # for ii in EC
    #     lawscc[ii]=μc[ii]#lawsccξ[ii]#
    # end
    return sample
end


@everywhere function sample1(μb,Σb,Ξ)
    S=Array{Float64}(undef,length(μb),length(Ξ))
    for ξ in Ξ
        S[:,ξ]=μb + randn(length(μb)).*Σb
    end
    return S
end

@everywhere function updatesamples2(Ξ,dhb,bb,μb,Σb,surface_extracted,surface_clusterblocks,lawsccξ,DH,blockcluster,EC,meanscenarioinsampleadp)
    μnew,Σnew = updatedistribution2(dhb,bb,μb,Σb,surface_extracted,surface_clusterblocks,lawsccξ,DH,blockcluster,EC)
    # dhb,bb,μb,Σb,surface_clusters,surface_clusterblocks,lawsccξ,DH,false)
    # = updatedistribution(μc,C,EC,ECτ,R,round.(lawsccξ,4))
    # Σ12iR22=Σ12invR22(C,EC,ECτ,R)
    if meanscenarioinsampleadp
        lawscc=μnew
    else
        lawscc=sample1(μnew,Σnew,Ξ)
    end
    return lawscc
end

function Elipsoid_rest(A,xi,xj,diagonal_matrix=true)
    c=0.
    ℓ=0
    bol=false
    while c>=1. || ℓ>length(xi)
        c=c + A[ℓ,ℓ] * (xi[ℓ]-xj[ℓ])^2
        ℓ=ℓ+1
    end
    if c>=1.
        bol=true
    elseif !diagonal_matrix && ℓ>length(xi)
        c= quad(A, xi-xj)
        if c>=1.
            bol=true
        end
    end
    return !bol
end

function krigingcovariance(Is,Js,Vs,𝓒,x,𝓐,A,bci)
    for i in bci
        k=1
        while k<=length(bci)
            j=bci[k]
            if Elipsoid_rest(A,x[i],x[j])
                push!(Is,i)
                push!(Js,j)
                push!(Vs,  𝓒 * exp(- norm(x[i]-x[j])/𝓐 ))
                if i!=j
                    push!(Is,j)
                    push!(Js,i)
                    push!(Vs,  𝓒 * exp(- norm(x[i]-x[j])/𝓐 ) )
                end
            end
            k=k+1
        end
    end
    mat=sparse(Is,Js,Vs)
    return mat#PDSparseMat(mat)#
end

#A=PDSparseMat(sparse(eye(3)/(9 * 𝓐^2)))
function kriging_on_blocks(bci,x,x_DH,DH,𝓐,A,dhb=[])#
    𝓒=var(DH)
    if dhb==[]
        for (i,pos) in enumerate(x)
            push!(dhb,Int32[])
            for (j,dh) in enumerate(x_DH)
                if Elipsoid_rest(A,pos,dh)
                    push!(dhb[i],j)
                end
            end
        end
    end
    Is=Int64[]
    Js=Int64[]
    Vs=Float64[]
    μB=Float32[]

    ΣDH=krigingcovariance(Is,Js,Vs,𝓒,x_DH,𝓐,A,1:length(x_DH))
    for (k,i) in enumerate(bci)
        Λdhc= PDMats.full(ΣDH[dhb[i],dhb[i]]) \ (𝓒 * exp.(
                [-norm(x[i]-x_DH[dh])/𝓐 for dh in dhb[i]]
            )
        )

        push!(μB,sum(Λdhc[𝐪]*DH[q] for (𝐪,q) in enumerate(dhb[i]) ))
    end
    # Is=Int64[]
    # Js=Int64[]
    # Vs=Float64[]
    return μB#,krigingcovariance(Is,Js,Vs,𝓒,x,𝓐,A,1:length(x))


end

function kriging_on_clusters(𝝹,clusterblocks,dhc,x,DH,x_DH,𝓐,x_C,wdir,numdhldr,Σbolean=[true,true])
    𝓒=var(DH)
    μC=zeros(𝝹)
    for (𝐜,c) in enumerate(clusterblocks)
        ΣP=𝓒*Matrix{Float64}(I, length(dhc[𝐜]), length(dhc[𝐜]))
        for 𝐩 =1:length(dhc[𝐜])
            p=x_DH[𝐩,:]
            for 𝐪=1:length(dhc[𝐜][(𝐩+1):end])
                q=x_DH[𝐪,:]
                ΣP[𝐩,𝐩+𝐪]=ΣP[𝐩+𝐪,𝐩]= round(𝓒 * exp(- norm(p-q)/𝓐 ),digits=6)
            end
        end
        # Σdhc=zeros(dhc[𝐜])
        # for t in c
        #     Σdhc=Σdhc + 𝓒 * exp.([-norm(x[t]-x_DH[dh])/𝓐 for dh in dhc[𝐜]]) / length(c)
        # end
        # Σdhc = try
        Σdhc =𝓒 * exp.([-norm(x_C[𝐜,:]-x_DH[dh,:])/𝓐 for dh in dhc[𝐜]]);
       # catch y
       #     println(y);
       #     println("size(x_C)=",size(x_C));
       #     println("size(x_DH)=",size(x_DH));
       #     println("𝐜=",𝐜);
       #     println("dhc[𝐜]=",dhc[𝐜]);
       #     println("x_DH[dhc[𝐜]]=",x_DH[dhc[𝐜]])
       #     exit();
       # end
        Λdhc=(Symmetric(ΣP\Matrix{Float64}(I, length(dhc[𝐜]) , length(dhc[𝐜]) )) )* Σdhc
        μC[𝐜] = μC[𝐜] + sum(Λdhc[𝐪]*DH[q] for (𝐪,q) in enumerate(dhc[𝐜]) )

    end
    # eye(𝝹)*𝓒
    ΣC=Matrix{Float64}(I, 𝝹, 𝝹)*𝓒
#    wwdir=replace(wdir,string(numdhldr)=>"")*"scenarios/sigma-"*string(Σbolean[1] ? 1 : 0)*"-"*string(Σbolean[2] ? 1 : 0)*".txt"
#     if !isfile(wwdir)
      for c1 in 1:length(clusterblocks)
          for c2 in (c1+1):length(clusterblocks)
              if Σbolean[1] && Σbolean[2]
                  ΣC[c1,c2]=ΣC[c2,c1]= 𝓒 * sum( exp(- norm(p-q)/𝓐 ) for p in clusterblocks[c1] for q in clusterblocks[c2])/ (length(clusterblocks[c1])*length(clusterblocks[c2]))
              elseif Σbolean[1]
                  ΣC[c1,c2]=ΣC[c2,c1]= 𝓒 * sum( exp(- norm(p-x_C[c2,:])/𝓐 ) for p in clusterblocks[c1])/ length(clusterblocks[c1])
              elseif Σbolean[2]
                  ΣC[c1,c2]=ΣC[c2,c1]= 𝓒 * sum( exp(- norm(x_C[c1,:]-q)/𝓐 ) for q in clusterblocks[c2])/ length(clusterblocks[c2])
              else
                  ΣC[c1,c2]=ΣC[c2,c1]= 𝓒 * exp(- norm(x_C[c1,:]-x_C[c2,:])/𝓐 )
              end
          end
      end
#          writedlm(wwdir,ΣC)
#     else
#         ΣC=readdlm(wwdir)
#     end
    return  μC,ΣC
end
#     μC = round.([-0.19802728404448947 0.24416056854364077 -1.0465822046709763 0.4213753736874735 -0.3800042735128818 -0.22291622257903235 -0.3533115553917617 -0.2994395496501188 -0.2516301953851409 0.007672679280185855 0.43754752876760733 -0.13763624200932734 -0.27935725373237535 -0.14367359139599525 -0.27198601885323453 -0.29378026878050534 0.13207708398582216 0.3540481987260351 -1.0475315730726442 -0.011794029986403465 -0.15493869729699855 -0.08028376644896082 -0.11769189149502764
#-0.17257597923896842 -0.043843301858089975 0.19814303119748616 0.3016844937434869 -0.15416786363822244 -0.2674590322353277 0.19579619715350952 -0.11834184266083883 -0.24602413404206253 0.26551106887902803 0.48643212703446653 -0.723057440459611 -0.34902238746272907 -0.0901384849758411 0.045972729756304266 -0.008016282675122766 -0.06552409443176667 0.24975602365207172 0.2219494784003843 0.2049050047644395 -0.3685477604227289 -0.14320445320215253 0.373853158461947 -0.014625979690987966
# -0.15078376827816883 ],round_val)

function sample2(μB,x,𝓒,𝓐,BB,Ξ)
    smpl=Array{Float64}(undef,length(μB),length(Ξ))
    sampled=Int32[]
    for b = 1:length(μB)
        smplI=intersect(sampled,BB[b])
        if length(smplI)==0 || b==1
            smpl[b,:]= μB[b]*ones(length(Ξ)) + 𝓒 * randn(length(Ξ))
        else
            # reduce(hcat,[ν[blockcluster[i],:] for i in bci ])'
            auxΣi= length(smplI)>1 ? inv(reduce(hcat,[Float64[ exp( -norm(x[b1]-x[b2])/𝓐 ) for b2 in smplI ] for b1 in smplI])) *
            Float64[ exp( -norm(x[b]-x[b1])/𝓐 ) for b1 in smplI ] : Float64[exp( -norm(x[b]-x[smplI[1]])/𝓐 )]
            for ξ in Ξ
                auxμ=μB[b] + auxΣi' * (smpl[smplI,ξ] - μB[smplI])
                auxσ=𝓒 - auxΣi' * Float64[ 𝓒 * exp( -norm(x[b]-x[b1])/𝓐 ) for b1 in smplI ]
                smpl[b,ξ] = auxμ + randn()*auxσ
            end
        end
        push!(sampled,b)
    end
    return smpl
end



#
# function kriging(𝝹,blockcluster,x,E_B,P,DH,x_DH,E_DH,dh_threshold,𝓐)
#     𝓒=var(DH)
#     μC=zeros(𝝹)
#     counterC=zeros(𝝹)
#     ΣP=𝓒*eye(length(x_DH))
#     for (𝐩,p) in enumerate(x_DH)
#         for (𝐪,q) in enumerate(x_DH[(𝐩+1):end])
#             ΣP[𝐩,𝐩+𝐪]=ΣP[𝐩+𝐪,𝐩]= 𝓒 * exp(- norm(p-q)/𝓐 )
#         end
#     end
#     e=1
#     while e<=dh_threshold
#         Σdht=Float32[]
#         dht=Int32[]
#         d=e
#         while d<=dh_threshold && E_B[e]==E_B[d]
#             push!(Σdht,𝓒 * exp(- norm(x[E_B[d]]-x_DH[E_DH[d]])/𝓐 ) )
#             push!(dht,E_DH[d])
#             d=d+1
#         end
#         Λb=(Symmetric(ΣP[dht,dht]\eye(length(dht))) )* Σdht
#         𝐜=blockcluster[E_B[e]]
#         μC[𝐜] = μC[𝐜] + sum(Λb[𝐪]*DH[q] for (𝐪,q) in enumerate(dht) )
#         counterC[𝐜]=counterC[𝐜]+1
#         e=d
#     end
#     for 𝐜 in 1:𝝹
#         μC[𝐜]=μC[𝐜]/counterC[𝐜]
#     end
#     counterC=zeros(𝝹,𝝹)
#     ΣC=eye(𝝹)*𝓒
#     while e <= length(E_B)
#         c1=blockcluster[E_B[e]]
#         c2=blockcluster[E_DH[e]]
#         ΣC[c1,c2] = ΣC[c1,c2] + 𝓒 * exp(- norm(x[E_B[e]]-x[E_DH[e]])/𝓐 )
#         counterC[c1,c2]=counterC[c1,c2]+1
#         e=e+1
#     end
#     for c1 in 1:𝝹
#         for c2 in (c1+1):𝝹
#             if ΣC[c1,c2]>0
#                 ΣC[c2,c1]=ΣC[c1,c2]=ΣC[c1,c2]/counterC[c1,c2]
#             elseif ΣC[c2,c1]>0
#                 ΣC[c1,c2]=ΣC[c2,c1]=ΣC[c2,c1]/counterC[c2,c1]
#             end
#         end
#     end
#     return  μC,ΣC
# end

# wdir="/home/tomas/Dropbox/DMP/codes/data/";
#
# NoEscenarios = 100;
#
#
# nblocks,edges,𝝹,bench,phase,wext,info,scenarioslaws,blockcluster,clusterblocks=initparamsDMP(wdir,true)
#
# μ,R=builddist(scenarioslaws,NoEscenarios,𝝹)
# Σ12iR22ξ2μ2,Σ12iR22=Σ12invR22ξ2μ2(μ,scenarioslaws[:,23],[117,114],R)
#
# upsampl=updatesamples(μ,transpose(scenarioslaws),[117,114],Σ12iR22ξ2μ2)
# μprim,Rprim=updatedistribution(μ,[117,114],R,Σ12iR22,Σ12iR22ξ2μ2)



# function builddistcluster(scen,NoEscenarios,𝝹)
#     μ = sum(scen[:,ξ] for ξ = 1:NoEscenarios)/NoEscenarios
#     for i =1:𝝹
#         for j = i:𝝹
#             push!(args,[ [ scen[i,k] for k=1:NoEscenarios ] ,  [ scen[j,k] for k=1:NoEscenarios ] ,μ[i],μ[j]])
#         end
#     end
#     a = SharedArray{Float32}(Int(𝝹*(𝝹+1)/2))
#     # @parallel
#     # a=zeros(Int(𝝹*(𝝹+1)/2))
#     for k = 1:Int(𝝹*(𝝹+1)/2)
#         a[k] = f(args[k])
#     end
#
#     k=1
#     R=zeros(Float32,𝝹,𝝹)
#     for i =1:𝝹
#         for j = i:𝝹
#             R[i,j]=R[j,i]=(a[k])
#             k=k+1
#         end
#     end
#     return μ,R
# end

#
# function Σ12invR22(C,EC,R)
#     sd=setdiff(C,EC)
#     invR22 = (R[EC,EC])\eye(EC)
#     return R[sd,EC] * invR22
# end
