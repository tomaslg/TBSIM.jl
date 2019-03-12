
#builds antecesors and predecesors
# input:the edges and number of clusters
function buildpredant(edges,𝝹)
    A=[]
    B=[]
    for i = 1:𝝹
        push!(A,Int64[])
        push!(B,Int64[])
    end
    for e =1:length(edges[1])
        push!(A[edges[1][e]],edges[2][e])
        push!(B[edges[2][e]],edges[1][e])
    end
    A,B
end

#Records different kind of structures: matrixes,multiarrays,values
function record(wdir,matrixes,multiarrays,values,name=["datamatrix","datamultiarray","datavalues"])
  for i = 1:length(matrixes)
    writedlm(string(wdir,"log/",name[1],"-",i,".txt"),matrixes[i])#,";")
  end
  for i = 1:length(multiarrays)
    open(string(wdir,"log/",name[2],"-",i,".txt"), "w") do f
        trs=""
        if length(multiarrays[i][1])>0
          trs=string(multiarrays[i][1][1])
          for j in multiarrays[i][1][2:end]
            trs=string(trs,";",j)
          end
        end
        for k in multiarrays[i][2:end]
          trs=string(trs,"\n")
          if length(k)>0
            trs=string(trs,k[1])
            for j in k[2:end]
              trs=string(trs,";",j)
            end
          end
        end
        write(f,trs)
     end
  end
  if length(values)>0 writedlm(string(wdir,"log/",name[3],".txt"),values) end
end



#reads different kind of structures: matrixes,multiarrays,values
function readrecord(wdir,nmatrixes,nmultiarrays,bvalues=true,names=["datamatrix","datamultiarray","datavalues"])
  matrixes=[]
  for i =1:nmatrixes
    # println(i)
    push!(matrixes,readdlm(string(wdir,"log/",names[1],"-",i,".txt"), (i<=2 ? Float64 : Int64) ))
  end
  multiarrays=[]
  for i=1:nmultiarrays
    push!(multiarrays,[])
    open(string(wdir,"log/",names[2],"-",i,".txt")) do f
        for (k,j) in enumerate(eachline(f))
          push!(multiarrays[i],Int64[])
          if j!=""
            io=split(j, r";");
            for ii in io
                push!(multiarrays[i][k],parse(Int,ii))
            end
          end
        end
    end
  end
  if bvalues
    values=readdlm(string(wdir,"log/",names[3],".txt"),Int64)
    for i = 1:length(values)
      values[i]=Int(values[i])
    end
    return matrixes,multiarrays,values
  else
    return matrixes,multiarrays
  end
end

#Removes the block that can not be reached in a feasible open-pit
function removeinfeasible(ξ,x,numlev,ldr,sideblocks=true)
  A=[0 0 1;0 1 1;1 0 1; 0 -1 1;-1 0 1]
  b=[numlev,800,800,0,0]
  y=[]
  ξnew=[Float64[] for nw=1:length(ξ[1,:])]
  k=0
  for k = 1:length(ξ[:,1])
    𝐱=x[k,:]
    bol=true
    for bl in [dot(A[i,:],𝐱) <= b[i] + ((sideblocks || i == 1) ? 0 : 20) for i = 1:length(b)]
      if !bl
        bol=false
        break
      end
    end
    if bol
      if sideblocks
        bol = ((𝐱[1]-5)/10)%(2^ldr)==0 && ((𝐱[2]-5)/10)%(2^ldr)==0
      else
        bol = ((𝐱[1]+20)/40)%(2^ldr)==0 && ((𝐱[2]+20)/40)%(2^ldr)==0
      end
    end
    if bol
      push!(y,𝐱)
      for nw=1:length(ξ[k,:])
        push!(ξnew[nw],ξ[k,nw])
      end
    end
  end
  return ξnew,y
end

#
function formclusters(x,includephaseprec,dim=[800, 800, 300],phases_length=[200,200,20])
  phases_number=[dim[i]/phases_length[i] for i=1:length(phases_length)]#[4,4,15]
  #build poliedral restrictions for the feasible blocks
  b=[ [ j*dim[k]/i for j=1:(i-1)] for (k,i) in enumerate(phases_number) ]
  #vcat([[1 0 0 ; 1 0 0 ; 1 0 0 ; 0 1 0 ; 0 1 0 ; 0 1 0]],
  A=[ [ k==1 ? Float16[1 0 0] : (k==2 ? Float16[0 1 0] : Float16[0 0 1])  for j=1:(i-1)] for (k,i) in enumerate(phases_number) ]
  # bench=[]

  nphase=[] #Phase id set for blocks
  blockcluster=Int16[] #each block's cluster id
  clusterblocks=[Int32[] for i=1:prod(phases_number)] #all blocks within each cluster
  for (j,𝐱) in enumerate(x)
    push!(nphase,ones(Int8,length(𝐱)))#for each phase(Easting,Northing,Elevation)
    for i=1:length(𝐱)
      k=1
      while k<=(phases_number[i]-1) && (A[i][k]*𝐱)[1] >= b[i][k]
        nphase[j][i] = nphase[j][i] + 1
        k=k+1
      end
    end
    cluster=Int(sum( (nphase[j][i]-1)*prod(phases_number[1:(i-1)]) for i=2:length(𝐱) ) + nphase[j][1])
    push!(clusterblocks[cluster],j)
    push!(blockcluster,cluster)
  end
  𝝹=0
  I=Int32[]
  for c in clusterblocks
    if length(c)>0
      𝝹=𝝹+1
      push!(I,𝝹)
    else
      push!(I,0)
    end
  end

  for j=1:prod(phases_number)
    i=Int32(prod(phases_number)-j+1)
    if I[i]==0
      splice!(clusterblocks,i)
    end
  end

  for i=1:length(blockcluster)
    blockcluster[i]=I[blockcluster[i]]
  end

  #At this point both mappings: from b \to c and from c \to b \in c are ready:
    # The former is blockcluster and the later clusterblocks.

  #As a result we form the edges that define cluster preferences
  edges=[Int32[],Int32[]]
  if includephaseprec[1] || includephaseprec[2] #Two types of precedences in the plane x,y
    for i=1:Int32(phases_number[2])
      for k=1:Int32(phases_number[1])
        cur=Int32((i-1)*phases_number[1] + k)
        if I[cur]>0# && I[cur - Int32(phases_number[1]*phases_number[2])]>0
          if includephaseprec[1] && k>1 && I[cur-1]>0
            push!(edges[1],cur - 1)
            push!(edges[2],cur )
          end
          if includephaseprec[2] && i>1 && I[cur-Int32(phases_number[1])]>0
            push!(edges[1],cur - Int32(phases_number[1]))
            push!(edges[2],cur )
          end
        end
      end
    end
  end
  for j=2:Int32(phases_number[end])
    for i=1:Int32(phases_number[2])
      for k=1:Int32(phases_number[1])
        cur=Int32((j-1)*phases_number[1]*phases_number[2] + (i-1)*phases_number[1] + k)
        if I[cur]>0# && I[cur - Int32(phases_number[1]*phases_number[2])]>0
          push!(edges[1],cur - Int32(phases_number[1]*phases_number[2]))
          push!(edges[2],cur)
          if includephaseprec[1] && k>1 && I[cur-1]>0
            push!(edges[1],cur - 1)
            push!(edges[2],cur )
          end
          if includephaseprec[2] && i>1 && I[cur-Int32(phases_number[1])]>0
            push!(edges[1],cur - Int32(phases_number[1]))
            push!(edges[2],cur )
          end
        end
      end
    end
  end

  #gets the pairs bech-phase that identifies each cluster
  benchphase=[[],[]]
  k=0
  for i=1:(Int32(phases_number[end]))
    for j=1:(Int32(phases_number[1]*phases_number[2]))
      k=k+1
      if I[k]>0
        push!(benchphase[1],i)
        push!(benchphase[2],j)
      end
    end
  end

  #These (x_C) are the positions of the clusters.
  # x_C=[copy(x[c[1]]) for c in clusterblocks]
  x_C=Array{Float64}(undef,length(clusterblocks),3);
  #The set of blocks that are on the surfaces of each cluster
  surface_clusterblocks=[]
  for c in 1:length(clusterblocks)
    x_C[c,end]=phases_length[end]*benchphase[1][c]
    x_C[c,2]=phases_length[2]*floor((benchphase[2][c]-1)/phases_number[1])
    x_C[c,1]=phases_length[1]*((benchphase[2][c]-1)%phases_number[1])
    push!(surface_clusterblocks,Int32[])
    y1=minimum([ x[i][1] for i in clusterblocks[c]])
    y2=maximum([ x[i][1] for i in clusterblocks[c]])
    w1=minimum([ x[i][2] for i in clusterblocks[c]])
    w2=maximum([ x[i][2] for i in clusterblocks[c]])
    z=maximum([ x[i][3] for i in clusterblocks[c]])
    for i in clusterblocks[c]
      if x[i][1]==y1 || x[i][1]==y2 || x[i][2]==w1 || x[i][2]==w2 || x[i][3]==z
        push!(surface_clusterblocks[c],i)
      end
    end
  end
  return benchphase,𝝹,clusterblocks,blockcluster,edges,x_C,surface_clusterblocks
end

#The following model is used to get the minimal set of index that define active constraints in the voronoi region
function get_redundant(A,b,i)
    m=Model(solver=mastersolver)
    @variable(m,y[h in 1:length(A[i,:])]>=0);
    @objective(m,:Min,
        -sum(y[h]*A[i,h] for h=1:length(A[i,:]))
        );
    @constraint(m,myc[k in [p for p in 1:length(b) if p!=i] ],
        sum(y[h]*A[k,h] for h=1:length(A[k,:]))  <= b[k]
    );
    solve(m,suppress_warnings=true);
    return (b[i]+getobjectivevalue(m)<0.)
end

#For each x, gets the x_DH (drill holes) indexes (store in dhb[x])
# such that define active constarints in the
# voronoi region given by x[i] at the center
# and x_DH arround. This function also returns (for each x)
# the block indexes that share  ρ×100% of the x_DH (drill holes) indexes
function intersect2(k9,l8,ρ)
  bol=true
  K0=0
  k9l8= union(k9,l8)
  L0=ρ*length(k9l8)
  # println("L0=",L0)
  for (𝔦,i) in enumerate(k9)
    for (𝔧,j) in enumerate(l8)
      if i==j
        K0=K0+1;
        if K0>=L0
          return true;
        end
      end
    end
  end
  if bol
    return false;
  end
end
function get_voronoi_indexes(x,x_DH)
    dhb=[]
    for i=1:length(x)
        ℓ1=ℓ2=1;
        bol=true;
        while ℓ2 <= length(x_DH)
            if bol && x_DH[ℓ2][3]>=x[i][3]-10
                bol=false;
                ℓ1=ℓ2
            elseif x_DH[ℓ2][3]>x[i][3]+10
                break
            end
            ℓ2=ℓ2+1
        end
        ℓ2=ℓ2-1
        ℓ1=ℓ1-1
        dist=[norm(x_DH[j+ℓ1]-x[i]) for j=1:(ℓ2-ℓ1)]
        s=sortperm(dist)
        neighbors=min(20,ℓ2-ℓ1)
        s=s[1:neighbors]
        A=Array{Float32}(undef,neighbors,length(x[i]))
        b=zeros(Float32,neighbors)
        for j=1:neighbors
            A[j,:]=x_DH[s[j]+ℓ1]-x[i]
            b[j]=(dot(x_DH[s[j]+ℓ1],x_DH[s[j]+ℓ1])-dot(x[i],x[i]))/2
        end
        push!(dhb,Int32[])
        let j=1
          while j<=length(b)
              if get_redundant(A,b,j)
                  push!(dhb[i],s[j]+ℓ1)
                  j=j+1
              else
                  A=A[[k for k=1:length(b) if k!=j],:]
                  b=b[[k for k=1:length(b) if k!=j]]
                  splice!(s,j)
              end
          end
        end
    end
    return dhb#=,bb=#
end

function get_voronoi_indexes_for_single_dh(x)
  dhb=[]; depth=Int32(1 + ((x[end][3]-5)/10) );
  for i=1:length(x)
    level=Int32(1 + ((x[i][3]-5)/10) );
    if (level>1)
      push!(dhb,Int32[level-1,level]);
    else
      push!(dhb,Int32[level])
    end
    if level<depth
      push!(dhb[i],level+1)
    end
  end
  return dhb#=,bb=#
end

function get_bb(x,ρ=30.,𝔨=3)
  bb=[];setk=[];cur_set=Int32[];prev=1;
  for i=1:length(x)
    level=Int32(1 + ((x[i][3]-5)/10) ); push!(bb,Int32[]);
    if level==prev
      push!(cur_set,i);
    else
      push!(setk,copy(cur_set));
      cur_set=Int32[]; prev=level;
    end
    for ℓ=(max(level-𝔨,1)):(level-1)
      for iprime in setk[ℓ]
        if norm(x[i]-x[iprime]) <= ρ
          push!(bb[i],iprime);
        end
      end
    end
  end
  return bb
end

#This is the constructor of the instance
function initparamsDMP(argnum,instance,wdir,includephaseprec,numlev=300,𝓐=50.,α=2.5,offlineupdate=true,ldr=0,num_drill_holes=400,voronoi_bool=true,xavier=false)
  name="Grilla.dat"
  pos=readdlm( replace(replace(replace(wdir,string(Int(num_drill_holes/(4^ldr)))=>""),instance=>""),replace(string(argnum),string(Int(num_drill_holes/(4^ldr)))=>"")=>"")*"scenarios/"*name )
  if xavier
    name="laws.scen"
    lws=readdlm( replace(replace(replace(wdir,string(Int(num_drill_holes/(4^ldr)))=>""),instance=>""),replace(string(argnum),string(Int(num_drill_holes/(4^ldr)))=>"")=>"")*"scenarios/"*name )
    DH=lws[1:12000]
    ξnew=lws[12001:end]
    name="simulaciones.scen"
    Ξnew=readdlm( replace(replace(replace(wdir,string(Int(num_drill_holes/(4^ldr)))=>""),instance=>""),replace(string(argnum),string(Int(num_drill_holes/(4^ldr)))=>"")=>"")*"scenarios/"*name )
  else
    lws=tbsim(2.5,pos,[],[],1);
    DH=lws[1:12000];
    ξnew=lws[12001:end];
    Ξnew=!voronoi_bool && !use_kriging ? α*log.(tbsim(2.5,pos[12001:end,:]),pos[1:12000,:],DH,100,false) : Array{Float32}(undef,0,0);
  end
  if  num_drill_holes==1
    bestdh= argmin([ sum([DH[i + j*400] for j=0:(Int(numlev/10 - 1)) ]) for i=1:400])#reduce(hcat,[ ν[blockcluster[i],:] for i in bci ])'
    selected_dh=Int32[bestdh]
  # elseif num_drill_holes==5
  #   Msel = zeros(20,20)
  #   Msel[1,1]=Msel[end,end]=Msel[end,1]=Msel[1,end]=Msel[10,10]=1
  #   selected_dh=Int32[]
  #   for i=1:20
  #     for j=1:20
  #       if Msel[i,j]>0
  #         push!(selected_dh,(i-1)*20 + j)
  #       end
  #     end
  #   end
  elseif num_drill_holes==5 || num_drill_holes==10 || num_drill_holes==20
    #if offlineupdate
    Msel = Matrix{Float64}(I, 20, 20);
    for i=1:20
      if i<=(20-num_drill_holes)/2. || i>(20+num_drill_holes)/2.#i%(num_drill_holes==5 ? 4 : ( num_drill_holes==10 ? 2 : 1) ) != 0
        Msel[i,i]=0
      end
    end
    Msel = Msel[shuffle(1:end), :]
    #writedlm(string(wdir,"log/Msel",num_drill_holes,".txt"),Msel)
    #else
    #  Msel=readdlm(string(wdir,"log/Msel",num_drill_holes,".txt"))
    #end
    selected_dh=Int32[]
    for i=1:20
      for j=1:20
        if Msel[i,j]>0
          push!(selected_dh,(i-1)*20 + j)
        end
      end
    end
  elseif num_drill_holes<400
    selected_dh=Int32[]
    available_dh=Int32[i for i=1:400]
    for i=1:num_drill_holes
      nwn=rand(available_dh)
      push!(selected_dh,nwn)
      splice!(available_dh,findfirst(i->i==nwn,available_dh))
    end
    sort!(selected_dh)
  else
    selected_dh=Int32[i for i in 1:(num_drill_holes/(4^ldr))]
  end
  selected_samples=Int32[]
  for j = 0:(Int(numlev/10) -1)
    for i in selected_dh
      push!(selected_samples,j*(num_drill_holes/(4^ldr)) + i)
    end
  end
  DH,x_DH=removeinfeasible(DH[selected_samples],pos[selected_samples#=1:12000=#,:],numlev,0,false)
  DH=DH[1]
  if xavier || (!voronoi_bool && !use_kriging)
    Ξnew,x=removeinfeasible(Ξnew,pos[12001:end,:],numlev,0)
  end
  ξnew,x=removeinfeasible(ξnew,pos[12001:end,:],numlev,0)
  ξnew=ξnew[1]
  # outputformat,datacoord,extreme_coord,ydata,simucoord,COORDS,RR=tbsim(hcat(x_DH...)',DH,hcat(x...)')#datacoord,ydata,simucoord)
  benchphase,𝝹,clusterblocks,blockcluster,edges,x_C,surface_clusterblocks=formclusters(x,includephaseprec,[800,800,numlev])
  # =formprecedences(𝝹,benchphase,includephaseprec)




  # E_DH1=0#Int32[]
  # E_DH2=0#Int32[]
  # E_B1=0#Int32[]
  # E_B2=0#Int32[]
  # E_DH,E_B=Int32[],Int32[]
  dh_threshold=0
  dhc,dhb=[],[]
  if true
    # TT = stdout # save original STDOUT stream
    # redirect_stdout()
    # redirect_stdout(TT) # restore STDOUT
    # for j = 1:length(ξnew)
    #   for i = 1:length(DH)
    #     h_DH_B= norm(x_DH[i] - x[j])
    #     if 0<h_DH_B && h_DH_B<=𝓐#3*
    #       dh_threshold=dh_threshold+1
    #       push!(E_B,j) # E_DH1=E_DH1+1 #
    #       push!(E_DH,i) # E_DH2=E_DH2+1 #
    #     end
    #   end
    # end
    #if voronoi_bool
    dhb#=,bb=#=num_drill_holes==1 ? get_voronoi_indexes_for_single_dh(x) : get_voronoi_indexes(x,x_DH)
    #end
    for c =1:𝝹
      push!(dhc,Int32[])
    end
    for e = 1:length(dhb)
      for ee in dhb[e]
        if !(ee in dhc[blockcluster[e]])
          push!(dhc[blockcluster[e]],ee)
        end
      end
    end
    # for j = 1:length(ξnew)
    #   for i = (j+1):length(ξnew)
    #     h_B_B= norm(x[i]-x[j])
    #     if (0<h_B_B && h_B_B<=3*𝓐
    #       && !(i in clusterblocks[blockcluster[j]]))
    #       # E_B1=E_B1+1
    #       # E_B2=E_B2+1
    #       push!(E_B,j)
    #       push!(E_DH,i)
    #     end
    #   end
    # end



    # BB=[]
    # if Σboolean
    #   for (i,pos1) in enumerate(x)
    #       push!(BB,Int32[])
    #       for (j,pos2) in enumerate(x)
    #           if Elipsoid_rest(PDSparseMat(sparse(eye(3)/(9 * 𝓐^2))),pos1,pos2)
    #               push!(BB[i],j)
    #           end
    #       end
    #   end
    # end
    # if voronoi_bool
   # record(wdir,[#=E_DH,E_B=#],[dhc,dhb#=,bb=#],[dh_threshold,ldr,numlev],["datamatrix","datamultiarray","datavalues"])
    # else
    #   record(wdir,[#=E_DH,E_B=#],[dhc],[dh_threshold,ldr,numlev],["datamatrix","datamultiarray","datavalues"])
    # end
  else
    RRd=readrecord(wdir,0,voronoi_bool ? 2 : 1,true,["datamatrix","datamultiarray","datavalues"])
    # E_DH,E_B=[Int32(𝐞) for 𝐞 in RRd[1][1]],[Int32(𝐞) for 𝐞 in RRd[1][2]]
    if voronoi_bool
      dhc,dhb#=,bb=#=RRd[2]
    else
      dhc=RRd[2][1]
    end
    dh_threshold,ldr_record,numlev_record=RRd[3]
    if ldr!=ldr_record || numlev!=numlev_record
      exit(2)
    end
    # 𝖊=1
    # while 𝖊<length(E_DH)
    #   if ((x[E_DH[𝖊][1]][1]-5)/10)%(2^ldr)==0 &&
    #      ((x[E_DH[𝖊][1]][2]-5)/10)%(2^ldr)==0
    #     𝖊=𝖊+1
    #   else
    #     deleteat!(E_DH, 𝖊)
    #   end
    # end
    # 𝖊= dh_threshold+1
    # while 𝖊<length(E_B)
    #   if ((x[E_B[𝖊]][1]-5)/10)%(2^ldr)==0 &&
    #      ((x[E_B[𝖊]][2]-5)/10)%(2^ldr)==0 &&
    #      ((x[E_DH[𝖊]][1]-5)/10)%(2^ldr)==0 &&
    #      ((x[E_DH[𝖊]][2]-5)/10)%(2^ldr)==0
    #     𝖊=𝖊+1
    #   else
    #     deleteat!(E_B, 𝖊)
    #     deleteat!(E_DH, 𝖊)
    #   end
    # end
  end
  Ω= ones(length(ξnew)) * 2650 * 192000/length(ξnew)
  Ωext= Float64[sum(Ω[j] for j in clusterblocks[i]) for i = 1:𝝹]
  pred,ant=buildpredant(edges,𝝹)
  bb=get_bb(x)


  return DH,x_DH,ξnew,Ξnew,x,benchphase,𝝹,
          clusterblocks,blockcluster,1:length(blockcluster),
          edges,Ω,Ωext,#=E_DH,E_B,=#dh_threshold,pred,ant,dhc,x_C,dhb,bb,surface_clusterblocks
  # return nblocks,bci,edges,𝝹,bench,phase,wext,infos[:,6],scenarioslaws,blockcluster,clusterblocks,numlev
end



# function get_surface(τ,𝝹,wx,k,phase,bench)
#     if τ>1
#         𝐭=(τ-1)
#         for i in 1:𝝹
#             j=𝝹+1-i;
#             if wx[j,𝐭]>0
#             #   println("Period=",𝐭,", extract bench=",bench[j],", phase=",phase[j]," \n")
#               if k[phase[j]]>-1 && bench[k[phase[j]]] > bench[j] k[phase[j]]=j end
#             end
#         end
#     end
#     return
# end

# function get_surface_excavated(τ,𝝹,wx,k,phase,bench,kk,pred,ant)
#         # get_surface(τ,𝝹,wx,k,phase,bench)
#         for i = 1:length(kk)
#             k[i]=kk[i]
#             while sum(wx[k[i],1:τ])>0
#                 #println("pred[k[i="*string(i)*"]="*string(k[i])*"]="*string(pred[k[i]]))
#                 BOL=true
#                 for j in pred[k[i]]
#                     if phase[j]==i
#                         k[i]=j;
#                         BOL=false
#                         break
#                     end
#                 end
#                 if BOL k[i]=-1; break; end
#             end
#             if k[i]!=-1 && i!=1 && sum(wx[ant[k[i]],1:τ])<length(ant[k[i]])
#                 k[i]=-1
#             end
#         end
#         return
#     end

#mean_and_std(scenarioslaws,2)
function get_surface_excavated(τ,𝝹,wx,k,edges)#,phase,bench,kk,pred,ant)
    # get_surface(τ,𝝹,wx,k,phase,bench)
    sumw=sum(wx[:,i] for i=1:τ)
    S=Int32[i for i=1:𝝹 if sumw[i]>0.99]
    e=1
    surface_extracted=Int32[]
    surface_unextracted=Int32[]
    while e<=length(edges[1]) && (edges[2][e] in S || !(edges[1][e] in S))  e=e+1 end
    while e<=length(edges[1])
        surface_extracted=union(surface_extracted,Int32[edges[1][e]])
        surface_unextracted=union(surface_unextracted,Int32[edges[2][e]])
        # push!(surface_unextracted,edges[2][e])
        if sum(sumw[ant[edges[2][e]]])>=length(ant[edges[2][e]])
            k[1 + (edges[2][e]-1)%length(k)]=edges[2][e]
        end
        e=e+1
        while e<=length(edges[1]) && #=(edges[1][e] in S || !(edges[2][e] in S))  &&=# (!(edges[1][e] in S) || edges[2][e] in S) e=e+1 end
    end
    if sum(sumw)<0.01
        k[1]=1
        push!(surface_unextracted,1)
    end
    # for i = 1:length(kk)
    #     k[i]=kk[i]
    #     while sum(wx[k[i],1:τ])>0.99
    #         #println("pred[k[i="*string(i)*"]="*string(k[i])*"]="*string(pred[k[i]]))
    #         BOL=true
    #         for j in pred[k[i]]
    #             if phase[j]==i
    #                 k[i]=j;
    #                 BOL=false
    #                 break
    #             else
    #                 k[i]=0
    #             end
    #         end
    #         if k[i]==0 || length(pred[k[i]]) == 0 break end
    #         # if BOL k[i]=-1; break; end
    #     end
    #     if k[i]>0 && i!=1 && sum(wx[ant[k[i]],1:τ])<length(ant[k[i]])
    #         k[i]=-1
    #     end
    # end
    return surface_extracted,S#,surface_unextracted
end

# function formprecedences(𝝹,benchphase,includephaseprec)
#   bench,phase=benchphase
#   edges= [Int32[],Int32[]];
#   for j =1:𝝹
#       if phase[j]==1 #The precedences in the first phase are only between benches
#         if bench[j]>bench[1]
#           for k =1:𝝹
#             if bench[k]+1 == bench[j] && phase[k]==phase[j]
#               push!(edges[1],k);
#               push!(edges[2],j);
#               break;
#             end
#           end
#         end
#       elseif bench[j]>bench[1]
#         #For any other phase the precedences are between benches and phases, except for the top level bench
#         for k =1:𝝹
#           if bench[k]+1 == bench[j] && phase[k]==phase[j]
#             push!(edges[1],k);
#             push!(edges[2],j);
#             break;
#           end
#         end
#         if includephaseprec
#           for k =1:𝝹
#             if bench[k] == bench[j] && phase[k]+1==phase[j]
#               push!(edges[1],k);
#               push!(edges[2],j);
#               break;
#             end
#           end
#         end
#       else #The top level bench only considers precedences between phases
#         if includephaseprec
#           for k =1:𝝹
#             if bench[k] == bench[j] && phase[k]+1==phase[j]
#               push!(edges[1],k);
#               push!(edges[2],j);
#               break;
#             end
#           end
#         end
#       end
#   end
#   return edges
# end
