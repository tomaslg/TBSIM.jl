# addprocs(3)
#
#
using Distributed,Distributions,JuMP
@everywhere using Gurobi##,
@everywhere env = Gurobi.Env()
@everywhere mastersolver = GurobiSolver(env,OutputFlag=0,Threads=1,PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0)#
@everywhere cbsolver = GurobiSolver(env,OutputFlag=0,Threads=1,MIPGap=0.,PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0)
include("initparamsDMP.jl")
include("builddist.jl")
# #
#
using TBSIM

function solveDEmodel(Captotalextraccion,Captotalproceso1,Ξ,T,τ,D,𝝹,r,Ωext,wp,scen,lawscb,Edges,clusterblocks,blockcluster,blocks)
  scenarioslaws=linearTprc(exp.(scen/α) /100,lawscb,wp,blockcluster,NoEscenarios,blocks)
  # linearTprc(scen,lawscb,wp,blockcluster,NoEscenarios,blocks)

  m = Model(solver = mastersolver);
  @variable(m, x[j=1:𝝹,t=1:(T-τ+1)] ,Bin);
  @variable(m, 0 <= y[i in blocks,t=1:(T-τ+1),d=1:D,ξ in Ξ] <= 1);
  for t=1:(T-τ+1)
    for c=1:𝝹
      for i in clusterblocks[c]
        for ξ in Ξ
          @constraint(m, sum(y[i,t,d,ξ] for d=1:D) == x[c,t]);
        end
      end
    end
  end

  for c=1:𝝹
    @constraint(m, sum(x[c,s] for s=1:(T-τ+1)) <= 1)
  end

  @objective(m, Max, sum( (1/(1+r))^(t+τ-2) * sum(y[i,t,d,ξ]*scenarioslaws[ i , D*ξ - (d%D) ] for d=1:D, i in blocks, ξ in Ξ)/length(Ξ) for t=1:(T-τ+1)))


#  print("Precedence constraints\n");
  @constraint(m,
    precedences[ee=1:length(Edges[1]), t=1:(T-τ+1)],
    sum(x[Edges[1][ee],s] for s=1:t) >= sum(x[Edges[2][ee],s] for s=1:t)
  );

#  print("KP constraints\n")
  @constraint(m,
    KP[t=1:(T-τ+1)],
    sum(x[c,t]*Ωext[c] for c=1:𝝹) <= Captotalextraccion
  );

  @constraint(m,
    KPr[t=1:(T-τ+1),ξ in Ξ],
    sum(y[i,t,2,ξ]*wp[i] for i in blocks) <= Captotalproceso1
  );

  status = solve(m)
  return getobjectivevalue(m),getvalue(x)
end




@everywhere function setxdict(xdict,x,𝝹,T)
  for ccc=1:𝝹
    for ttt=1:T
      xdict[(ccc,ttt)]=x[ccc,ttt];
    end
  end
end

@everywhere function definemaster(ʖ,ℓ,T,τ,Edges,Ωext,Captotalextraccion,xdict,freshstart,lowerbound=0.,cbsolv=false)
  if cbsolv
    m = Model(solver = cbsolver);
  else
    m = Model(solver = mastersolver);
  end
  #@variable(m, x[i=1:10], start=(i/2)
  # if string(mastersolver)[1:6] == "Gurobi"
  if freshstart
    @variable(m, x[j=1:length(ʖ),t=1:(T-τ+1)] ,Bin);
  else
    @variable(m, x[j=1:length(ʖ),t=1:(T-τ+1)] ,Bin, start=Int(round( min(xdict[(j,t+τ-1)],ℓ[j]*1.0) )));
  end
  @variable(m, lowerbound <=z<= sum(Ωext) * 10^5)
  @objective(m,:Max,z);
  @constraint(m,
    ext_once_every_cluster[i=1:length(ʖ)] ,
    sum(x[i,ϕ] for ϕ=1:(T-τ+1)) <= 1
  );

  # print("Precedence constraints\n");
  @constraint(m,
    precedences[ee=1:length(Edges[1]), t=1:(T-τ+1)],
    sum(x[ℓ[Edges[1][ee]],s] for s=1:t) >= sum(x[ℓ[Edges[2][ee]],s] for s=1:t)
  );

  # print("KP constraints\n")
  @constraint(m,
    KP[t=1:(T-τ+1)],
    sum(x[i,t]*Ωext[ʖ[i]] for i=1:length(ʖ)) <= Captotalextraccion
  );

  return m,x,z,ext_once_every_cluster,precedences,KP

end

#Returns a index list following decreasing order wrt each scenario law
@everywhere function buildSP(scenarioslaws,NoEscenarios,D,bci,Ω)
    SP=[]
    for ξ =1:NoEscenarios
        Α=Float64[]
        for α = 1:length(scenarioslaws[:,ξ*D ])
            push!(Α,(scenarioslaws[α,ξ*D - 1]-scenarioslaws[α,ξ*D ])/Ω[α])
        end
        sp=sortperm(Α)
        push!(SP,intersect(sp,bci))
    end
    SP
end


#Calculates the value of the solution of the CKP for all scenarios
#If dual true, also returns dual multipliers
@everywhere function solveCKP(w,ωp,scenarioslaws,SP,D,Captotalproceso1,blockcluster,dual=false,ʖ=Int64[],ℓ=Int64[])#,Trpefr𝛰,Trξ)
    Qx=Float64(0)
    if length(ʖ)==0
      ʖ=ℓ=union(blockcluster)
    end
    if dual
        λ=zeros(length(blockcluster),1)
        μ=zeros(1)
        threshold=Int64[]
        meet_cap_constraint=Bool[]
    end
    for ξ =1:1
        sp=SP
        b=length(sp)
        while b>=1
            if !(blockcluster[sp[b]] in ʖ)
                splice!(sp,b)
            end
            b=b-1
        end
        remcapacity = Captotalproceso1[1]
        profit=Float64(0)
        id=1
        if dual
            push!(threshold,0)
            push!(meet_cap_constraint,false)
        end
        while id <= length(sp) && remcapacity - ωp[ sp[id] ]*w[ℓ[blockcluster[sp[id]]]] > 10. ^(-4.)
            if w[ℓ[blockcluster[sp[id]]]] > 10. ^(-4.)
                # println("scenarioslaws[ sp[id=",id,"] , ξ*D ]*w[blockcluster[ sp[id]=",sp[id]," ]] >=
                #     scenarioslaws[ sp[id] , (ξ*D -1) ]*w[ℓ[blockcluster[sp[id]]]]=",
                #     scenarioslaws[ sp[id] , ξ*D ]*w[ℓ[blockcluster[sp[id]]]],">=",
                #     scenarioslaws[ sp[id] , (ξ*D -1) ]*w[ℓ[blockcluster[sp[id]]]])
                if scenarioslaws[ sp[id] , ξ*D ] >= scenarioslaws[ sp[id] , (ξ*D -1) ]
                    profit=profit + scenarioslaws[ sp[id] , ξ*D ]*w[ℓ[blockcluster[sp[id]]]]
                else
                    break;
                end
                remcapacity=remcapacity - ωp[ sp[id] ]*w[ℓ[blockcluster[sp[id]]]]
            end
            id=id+1
        end
        # println(id <= length(sp) )
        # println(scenarioslaws[ sp[id] , ξ*D ]*w[ℓ[blockcluster[sp[id]]]] >= scenarioslaws[ sp[id] , (ξ*D-1) ]*w[ℓ[blockcluster[sp[id]]]])
        # println((remcapacity - ωp[ sp[id] ]*w[ℓ[blockcluster[sp[id]]]]<0))
        if dual
          if id <= length(sp)
            threshold[ξ]=id
          else
            threshold[ξ]=id-1
          end
        end
        if id <= length(sp) &&
                scenarioslaws[ sp[id] , ξ*D ]*w[ℓ[blockcluster[sp[id]]]] >= scenarioslaws[ sp[id] , (ξ*D-1) ]*w[ℓ[blockcluster[sp[id]]]] &&
                    (remcapacity - ωp[ sp[id] ]*w[ℓ[blockcluster[sp[id]]]]<0)
            if dual
                meet_cap_constraint[ξ]=true
            end
            profit=profit + scenarioslaws[ sp[id] , ξ*D ]*w[ℓ[blockcluster[sp[id]]]]* remcapacity/ωp[ sp[id] ]
            profit=profit + scenarioslaws[ sp[id] , ξ*D - 1]*w[ℓ[blockcluster[sp[id]]]] * (1 - remcapacity/ωp[ sp[id] ])
            id=id+1
        end
        while id <= length(sp)
            if w[ℓ[blockcluster[sp[id]]]] > 10. ^(-4.)
                profit=profit + scenarioslaws[ sp[id] , ξ*D - 1 ]
            end
            id=id+1
        end
        #if ξ == Trξ push!(Trpefr𝛰,profit) end
        Qx=Qx + profit
        # println("threshold=",threshold)
        if dual
          # thresh=findfirst(threshold[ξ],sp)
            for i=1:length(sp)
                di=length(sp)+1-i
                if di<threshold[ξ]
                    λ[sp[di],ξ]=scenarioslaws[sp[di],ξ*D ] - ωp[ sp[di] ]*μ[ξ]
                elseif di==threshold[ξ]
                  λ[sp[di],ξ]= scenarioslaws[sp[di],ξ*D - 1]
                    if meet_cap_constraint[ξ]
                        μ[ξ]=(#scenarioslaws[sp[di],ξ*D] - scenarioslaws[sp[di],ξ*D -1]
                              #~it fails because there must be a space after the - before the 1
                              #   syntax: missing separator in array expression
                              scenarioslaws[sp[di],ξ*D] - λ[sp[di],ξ]
                            )/ωp[ sp[di] ]
                    end
                else
                    λ[sp[di],ξ]= scenarioslaws[sp[di],ξ*D - 1]
                end
                # if abs((λ[sp[di]]-λbenm[sp[di]])*w[blockcluster[ sp[di] ]])>10. ^(-4.)
                #   println("λ[sp[di=",di,"]]-λbenm[sp[di]=",sp[di],"]=",(λ[sp[di]]-λbenm[sp[di]])*w[blockcluster[ sp[di] ]])
                # end
            end
            # if abs(μ[ξ]-μbenm[ξ])>10. ^(-4.)
            #   println("μ-μbenm",μ[ξ]-μbenm[ξ])
            # end
        end
    end
    if dual
        return Qx,λ,μ
    else
        return Qx
    end
end




@everywhere function solveDMPbend(arg,xdict=Dict{Tuple{Int32,Int32},Float64}(),freshstart=false,givesolution=false,parallel=false,cbsolv=false)
  Captotalextraccion,Captotalproceso1,#= Captotalproceso1=Captotalproceso1[1]=#
  Ωext,ν,Ω,edges,blockcluster,bci,clusterblocks,wx,Ξ,τ,#=τ=1,=#T,D,r,angulocut,U,α=arg
  NoEscenarios=length(Ξ)
  scen=ν#Array{Float64}(length(bci),NoEscenarios)
  # for i in bci
  #     scen[i,:]=ν[i,:]
  # end
  # println("length(scen[:,1])=",length(scen[:,1]))
  # println("2*NoEscenarios=",2*NoEscenarios)
  scenarioslaws=linearTprc(exp.(scen/α) /100,Ω,blockcluster,NoEscenarios,bci)
  sp=buildSP(scenarioslaws,NoEscenarios,D,bci,Ω)
  𝝹=length(clusterblocks)

  Edges = [Int64[],Int64[]]
  if τ==1
    Edges=edges
  else
    for e in 1:length(edges[1])
      bol=true
      if sum(wx[edges[1][e],ϕ] for ϕ = 1:(τ-1)) > 0 || sum(wx[edges[2][e],ϕ] for ϕ = 1:(τ-1)) > 0
        bol=false
      end
      if bol
        push!(Edges[1],edges[1][e])
        push!(Edges[2],edges[2][e])
      end
    end
  end
#  if 0 in Edges[1] || 0 in Edges[2] print("Edges=",Edges) end
  ʖ = Int64[]
  ℓ = Int64[]
  blocks=Int32[]
  let
    j=1
      for i=1:𝝹
        if τ==1 || sum(wx[i,1:(τ-1)])<0.01
          push!(ʖ,i)
          push!(ℓ,j)
          blocks=union(blocks,clusterblocks[i])
          j=j+1
        else
          push!(ℓ,0)
        end
    end
  end


  if length(ʖ)==0
    if givesolution
      return zeros(𝝹,T-τ+1)
    else
      return 0.0
    end
  end

  blocks = setdiff(bci,setdiff(bci,blocks))

  # for c = 1:length(ʖ)
  #   for b in clusterblocks[ʖ[c]]
  #     push!(blocks,b)
  #   end
  # end

  # if angulocut
  #   U = solvecristalballmodel(Captotalextraccion,Captotalproceso1,NoEscenarios,T,τ,D,ʖ,ℓ,r,wext,info[:,6],scenarioslaws,Edges,clusterblocks,blocks)
  # end
  owe1 = 0.;
  if !freshstart
    for ττ=τ:T
      for ξ=1:NoEscenarios
         owe1=owe1 + solveCKP([ min(xdict[(j,ττ)],ℓ[j]*1.0)  for j=1:𝝹],Ω,scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^(ττ) /length(Ξ),sp[ξ],D,Captotalproceso1,blockcluster)
  #    sub,λ,μ=solveCKP(solτ[t+1],Ω,scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^t /length(Ξ),SP,
  #                            1,D,Captotalproceso1*(T-τ),blockcluster,true,ʖ,ℓ)
      end
    end
  end
  m,x,z,ext_once_every_cluster,precedences,KP=definemaster(ʖ,ℓ,T,τ,Edges,Ωext,Captotalextraccion,xdict,freshstart,owe1,cbsolv);
  # solτ=[]
  # for t = 0:(T-τ)
  #   if ℓ[j]>0 push!(solτ,ones(length(ʖ))) end
  # end
  # firstcut=true
  #U=sum(abs(scenarioslaws[j,ξ]) for j=1:length(scenarioslaws[:,1]), ξ=1:length(scenarioslaws[1,:]) )


  function cut(cb)
    # println("----\nInside cut callback")
    masterval=getvalue(z);
    sol=getvalue(x)
    if angulocut aexpr=zero(AffExpr) end
    bexpr=zero(AffExpr)
    solval=0;
    # if !firstcut
      solτ=[]
    # end
    cf=zeros(length(ʖ),T-τ+1)
    if parallel
      slaves=Future[]
      for t = 0:(T-τ)
        # if !firstcut
        push!(solτ,Float64[])
        for j in ʖ
          if ℓ[j]>0 push!(solτ[t+1],sol[ℓ[j],t+1]) end
        end
        # elseif t==T-τ
        #   firstcut=false
        # end


        for ϕ=1:length(Ξ)
          ξ=Ξ[ϕ]
          # sub,λ,μ = definesubproblem(solτ[t+1], scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^t /length(Ξ) ,blockcluster,blocks,ℓ,Ω,Captotalproceso1)
          # solval=solval + sub;# * (1+r)^t );
          # solveCKP(w,ωp,scenarioslaws,SP,NoEscenarios,D,Captotalproceso1,blockcluster,dual=false)
          # (solτ[t+1],Ω,scenarioslaws,SP,NoEscenarios,D,Captotalproceso1[1],blockcluster)* δ^(τ-1) )
          # solval1
          push!(slaves,@spawn solveCKP(solτ[t+1],Ω,scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^t /length(Ξ),sp[ξ],
                            D,Captotalproceso1,blockcluster,true,ʖ,ℓ))
        end
      end
      result=[]#fetch(𝐳)
      while length(slaves)>0
          push!(result,fetch(pop!(slaves)))
      end
      kres=0
      for t = 0:(T-τ)
        for ϕ=1:length(Ξ)
          kres=kres+1
          sub,λ,μ=result[kres]
          solval=solval + sub
          # if abs(sub-
          #   solval1[1])>(10. ^(-4.))
          #   println("procedure is not optimal to deliver: ",sub-solval1[1]- solval1[3][1]*Captotalproceso1)
          # #   println("solval1= ",solval1)
          # #   println("μvsμ=",μ - solval1[3][1] )
          # #   println("λvsλ=",Float64[λ[o] for o=1:length(λ)]-solval1[2][:,1])
          # end
          bexpr +=  μ[1]*Captotalproceso1# + sum(solτ[t+1][ℓ[blockcluster[b]]] for b in blocks) )
          for c =1:length(ʖ)
            for b in clusterblocks[ʖ[c]]
              cf[c,t+1]+=λ[b]
            end
          end
        end

        for c =1:length(ʖ)
          push!(bexpr,cf[c,t+1],x[c,t+1])
        end

        if angulocut
          for c =1:length(ʖ)
            #aexpr = aexpr + (solτ[t+1][c]==1 ? -1 : 1)*x[c,t+1]
            push!(aexpr, solτ[t+1][c]==1 ? -1.0 : 1.0, x[c,t+1])
          end
          aexpr= aexpr+ sum(solτ[t+1])
        end
      end
    else
      for t = 0:(T-τ)
        # if !firstcut
          push!(solτ,Float64[])
          for j in ʖ
            if ℓ[j]>0 push!(solτ[t+1],sol[ℓ[j],t+1]) end
          end
        # elseif t==T-τ
        #   firstcut=false
        # end


        for ϕ=1:length(Ξ)
          ξ=Ξ[ϕ]
          # sub,λ,μ = definesubproblem(solτ[t+1], scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^t /length(Ξ) ,blockcluster,blocks,ℓ,Ω,Captotalproceso1)
          # solval=solval + sub;# * (1+r)^t );
          # solveCKP(w,ωp,scenarioslaws,SP,NoEscenarios,D,Captotalproceso1,blockcluster,dual=false)
          # (solτ[t+1],Ω,scenarioslaws,SP,NoEscenarios,D,Captotalproceso1[1],blockcluster)* δ^(τ-1) )
          # solval1
          sub,λ,μ=solveCKP(solτ[t+1],Ω,scenarioslaws[ : , (D*ξ - 1):(D*ξ) ] * (1/(1+r))^t /length(Ξ),sp[ξ],
                            D,Captotalproceso1,blockcluster,true,ʖ,ℓ)
          solval=solval + sub
          # if abs(sub-
          #   solval1[1])>(10. ^(-4.))
          #   println("procedure is not optimal to deliver: ",sub-solval1[1]- solval1[3][1]*Captotalproceso1)
          # #   println("solval1= ",solval1)
          # #   println("μvsμ=",μ - solval1[3][1] )
          # #   println("λvsλ=",Float64[λ[o] for o=1:length(λ)]-solval1[2][:,1])
          # end
          bexpr +=  μ[1]*Captotalproceso1# + sum(solτ[t+1][ℓ[blockcluster[b]]] for b in blocks) )
          for c =1:length(ʖ)
            for b in clusterblocks[ʖ[c]]
              cf[c,t+1]+=λ[b]
            end
          end
        end

        for c =1:length(ʖ)
          push!(bexpr,cf[c,t+1],x[c,t+1])
        end

        if angulocut
          for c =1:length(ʖ)
            #aexpr = aexpr + (solτ[t+1][c]==1 ? -1 : 1)*x[c,t+1]
            push!(aexpr, solτ[(t+1)][c]==1 ? -1.0 : 1.0, x[c,t+1])
          end
          aexpr= aexpr+ sum(solτ[t+1])
        end
      end
    end
    # print("solval= \n =")
    # println(solval)
    # print("masterval= \n =")
    # println(masterval)
    if masterval-0.0001<=solval<=masterval+0.0001
      return
    end
    # println("Adding sample average scenarios elimination cut")
    # println("solval=",solval)
    if solval!=-Inf
      if angulocut
        aexpr= aexpr*(U-solval) + solval
        push!(aexpr,-1.0, z)
        @lazyconstraint(cb,aexpr >= 0)#@usercut(cb, aexpr >= 0)#
      end
      #bexpr=bexpr-z
      push!(bexpr,-1.0,z)
    end

    @lazyconstraint(cb, bexpr >= 0)#@usercut(cb, bexpr >= 0)#
    # push!(bexpr,-1.0,z); @constraint(m, bexpr >= 0); solve(m);
  end

  addlazycallback(m,cut)#,fractional=true)#,localcut=true)#addcutcallback(m, cut)


  solve(m)
  gtvx=getvalue(x)
  gtov=getobjectivevalue(m)

  if freshstart
    setxdict(xdict,gtvx,𝝹,T);
  end
  if givesolution
    return gtvx
  else
    return gtov#* (1/(1+r))^(τ-1)
  end

end

arg=[Captotalextraccion,Captotalproceso1,Ωext,ν,Ω,edges,blockcluster,bci,clusterblocks,wx,Ξ,τ,#=τ=1,=#T,D,r,angulocut,U,α]
xdict=Dict{Tuple{Int32,Int32},Float64}()
freshstart=false
givesolution=false
parallel=false
cbsolv=false
solveDMPbend(arg,xdict,freshstart,givesolution,parallel,cbsolv)

# function definesubproblem(sol,ω,blockcluster,blocks,ℓ,weights,Captotalproceso1)
#   sb=Model(solver=subproblemsolver);
#   @variable(sb,μ>=0);
#   @variable(sb,λ[b in blocks]);
#   @objective(sb,:Min,
#     Captotalproceso1*μ +
#     sum(λ[b]* (max(sol[ℓ[blockcluster[b]]],0.)) for b in blocks )
#   );
#   @constraint(sb,
#     diet1[b in blocks],
#     λ[b] + weights[b]*μ  >= ω[b,2]
#   );
#   @constraint(sb,
#     diet2[b in blocks],
#     λ[b] >= ω[b,1]
#   );
#   sta=solve(sb);
#   if sta==:Unbounded
#     print(sb)
#     writeLP(sb,string(wdir,"out"))
#   end
#   return getobjectivevalue(sb),getvalue(λ),getvalue(μ)
# end

# @everywhere function definesubproblem(sol,ω,clusterblocks,blocks,ʖ,weights,Captotalproceso1,D)
#
#   sb=Model(solver=subproblemsolver);
#   @variable(sb,λp>=0);
#   @variable(sb,λu[b in blocks,d=1:D] >= 0);
#   @variable(sb,υ[b in blocks]);
#   @objective(sb,:Min,
#     Captotalproceso1*λp +
#     sum(λu[b,d] for b in blocks,d=1:D ) +
#     sum( sum(υ[b] for b in clusterblocks[ʖ[c]]) * sol[c] for c =1:length(ʖ))
#   );
#   @constraint(sb,
#     diet1[b in blocks],#[i for c in ʖ, i in clusterblocks[c]]],
#     υ[b] + λu[b,1] + weights[b]*λp >= ω[b,1]
#   );
#   @constraint(sb,
#     diet2[b in blocks],
#     υ[b] + λu[b,2]  >= ω[b,2]
#   );
#   solve(sb);
#   return getobjectivevalue(sb),getvalue(λp),getvalue(λu),getvalue(υ)
# end
