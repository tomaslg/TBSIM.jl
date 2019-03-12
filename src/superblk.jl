#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function superblk(datacoord,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz);

#----------------------------
# Set up super-block strategy
#----------------------------

# Default parameters

    MAXSBX = 21;
    MAXSBY = 21;
    MAXSBZ = 21;
    ndata = size(datacoord,1);

    # Establish the number and size of the super blocks

    nxsup   = min(MAXSBX,nx);
    nysup   = min(MAXSBY,ny);
    nzsup   = min(MAXSBZ,nz);
    xsizsup = nx*xsiz/nxsup;
    ysizsup = ny*ysiz/nysup;
    zsizsup = nz*zsiz/nzsup;
    xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup;
    ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup;
    zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup;

    # Assign the data to a super block

    ix = max.(1,min.(nxsup,Int.(floor.((datacoord[:,1] .- xmnsup)/xsizsup .+ 1.5))));
    iy = max.(1,min.(nysup,floor.((datacoord[:,2] .- ymnsup)/ysizsup .+ 1.5)));
    iz = max.(1,min.(nzsup,floor.((datacoord[:,3] .- zmnsup)/zsizsup .+ 1.5)));
    index = Int.(ix + (iy .- 1)*nxsup + (iz .- 1)*nxsup*nysup);
    # println("index=",index)

    #  Accumulate how many data belong to each super block

    nisb = Int.(zeros(nxsup*nysup*nzsup));
    for i = 1:ndata
      nisb[Int64(index[i])] = nisb[Int64(index[i])] + 1;
    end
    nisb = vcat(0,cumsum(nisb));
    # println("nisb=",nisb)

    # Sort by ascending super block number

    # [tmp,I] = sort(index);
    I=sortperm(index);
    return I,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup
end
