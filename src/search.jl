#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function Search(datacoord,datavalue,coord,search_rotationmatrix,octant,ndata,nxsup,nysup,
                                 nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb)

#-------------------------------------------------------------------------------------------
# Search for the data located near location "coord" according to the neighborhood parameters
# (search_rotationmatrix,octant,ndata), using the superblock search strategy
#-------------------------------------------------------------------------------------------


  if isempty(datacoord)
    return [],[];
  end


# First selection of data according to super-block search
#--------------------------------------------------------

# Indices of the super-block containing the node to simulate

  ix = max(1,min(nxsup,floor((coord[1]-xmnsup)./xsizsup + 1.5)));
  # println("ix=",ix);
  iy = max(1,min(nysup,floor((coord[2]-ymnsup)./ysizsup + 1.5)));
  # println("iy=",iy);
  iz = max(1,min(nzsup,floor((coord[3]-zmnsup)./zsizsup + 1.5)));
  # println("iz=",iz);

# Indices of the acceptable super-blocks for searching the data

  ix =  ix .+ ixsbtosr;#length(ixsbtosr)>0 ? () : ix;
  # println("ix=",ix);
  iy = iy .+ iysbtosr;#length(iysbtosr)>0 ? () : iy;
  # println("iy=",iy);
  iz = iz .+ izsbtosr;#length(izsbtosr)>0 ? () : iz;
  # println("iz=",iz);
  I = findall( (ix.>=1).&(iy.>=1).&(iz.>=1).&(ix.<=nxsup).&(iy.<=nysup).&(iz.<=nzsup) );
  # println("I=",I);
  ind_block = Int.(ix[I] + (iy[I] .- 1)*nxsup + (iz[I] .- 1)*nxsup*nysup);
  # println("ind_block=",ind_block);

# Indices of acceptable data

  n_block = length(ind_block);
  # println("n_block=",n_block)
  st = Int.(nisb[ind_block] .+ 1);
  # println("st=",st)
  finish = Int.(nisb[ind_block .+ 1]);
  # println("finish=",finish)
  index = Int.(vcat([0],cumsum(finish - st .+ 1)));
  # println("index=",index)
  I = zeros(Int32,index[n_block+1]);
  for i = 1:n_block
    # try
    I[(index[i]+1):index[i+1]] = (st[i]):(finish[i]);
    # catch ee
    #   println("I[(index[i]+1):index[i+1]] = (st[i]):(finish[i])");
    #   println("length(I)=",length(I));
    #   println("i=",i);
    #   println("index[i]+1=",index[i]+1);
    #   println("index[i+1]=",index[i+1]);
    #   println("st[i]=",st[i]);
    #   println("finish[i]=",finish[i]);
    #   println(I[(index[i]+1):index[i+1]]," = ",(st[i]):(finish[i]));
    #   println(ee);
    #   exit(1);
    #   # if isa(ee, DomainError)
    #   #     sqrt(complex(x[2], 0))
    #   # elseif isa(y, BoundsError)
    #   #     sqrt(x)
    #   # end
    # end
          # if n_block>0
          #
          #     if length(I)>0
          #       if index[i]+1<=length(I)
          #       end
          #       # try
          #       # catch ee
          #       #    if isa(ee,BoundsError)
          #       #      println("index[i]+1:index[i+1]=",index[i]+1:index[i+1]);
          #       #      println("(st[i]):(finish[i])=",(st[i]):(finish[i]));
          #       #      println("i=",i);
          #       #      println("size(I)=",size(I));
          #       #      println("finish=",finish);
          #       #      println("st=",st);
          #       #      println("cumsum(finish - st + 1)=",cumsum(finish - st + 1));
          #       #      println("n_block=",n_block);
          #       #
          #       #    end
          #       #    error(ee);
          #       # end
          #     end
          #   end
          # else
          #   I=[];
  end

# Select the acceptable neighboring data
  # I=Int.(I);
  datacoord_i = datacoord[I,:];
  residuals_i = datavalue[I,:];


# Second selection of data according to radius, angles and octant
#----------------------------------------------------------------

  n = size(datacoord_i,1);
  if (n == 0)
    return [],[];
  end
  # println("size(I)=",size(I))
  # println("size(datacoord_i)=",size(datacoord_i))
# Compute the reduced distances between data and location to estimate
  # println("datacoord_i - ones(n,1)*coord'=",datacoord_i," - ",ones(n,1)*coord')
  deltacoord = datacoord_i - ones(n,1)*coord;
  deltacoord = (deltacoord*search_rotationmatrix)';

# Create flags to indicate the belonging to angular sectors

  if (octant == 0)
    nsector = 1;
    flag=Array(nsector,n);
    flag[1,:] = [1:n];
  else
    nsector = 8;
    flag=Array{Int16}(undef,nsector,length(deltacoord[1,:]));
    flag[1,:] = (deltacoord[1,:] .> 0) .& (deltacoord[2,:] .> 0) .& (deltacoord[3,:] .> 0);
    flag[2,:] = (deltacoord[1,:] .> 0) .& (deltacoord[2,:] .> 0) .& (deltacoord[3,:] .<= 0);
    flag[3,:] = (deltacoord[1,:] .> 0) .& (deltacoord[2,:] .<= 0) .& (deltacoord[3,:] .> 0);
    flag[4,:] = (deltacoord[1,:] .> 0) .& (deltacoord[2,:] .<= 0) .& (deltacoord[3,:] .<= 0);
    flag[5,:] = (deltacoord[1,:] .<= 0) .& (deltacoord[2,:] .> 0) .& (deltacoord[3,:] .> 0);
    flag[6,:] = (deltacoord[1,:] .<= 0) .& (deltacoord[2,:] .> 0) .& (deltacoord[3,:] .<= 0);
    flag[7,:] = (deltacoord[1,:] .<= 0) .& (deltacoord[2,:] .<= 0) .& (deltacoord[3,:] .> 0);
    flag[8,:] = (deltacoord[1,:] .<= 0) .& (deltacoord[2,:] .<= 0) .& (deltacoord[3,:] .<= 0);

  end

# Select the neighboring data

  I = zeros(Int32,ndata*nsector,1);
  k = 0;

  for i = 1:nsector

    index = findall( xy->(xy > 0), flag[i,:]);#find(flag(i,:) > 0);

    if ~isempty(index)

      squareddistp = sum(deltacoord[:,index].^2,dims=1)';
      squareddist = Float64[];
      for ii in squareddistp push!(squareddist,ii) end
      # Sort the data by increasing distance
      J = sortperm(squareddist);
      sorteddist=squareddist[J];
      index = index[J];

      # Discard the data located beyond the radius
      J = findall( xy->(xy < 1), sorteddist);#find(sorteddist<1);
      n = min(ndata,length(J));
      index = index[1:n];
      I[k+1:k+n] = index;
      k = k+n;

    end

  end

  I = I[1:k];
  datacoord_i = datacoord_i[I,:];
  residuals_i = residuals_i[I,:];

  return datacoord_i,residuals_i
end
