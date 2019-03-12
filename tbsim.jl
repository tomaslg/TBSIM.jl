module TBSIM

  export tableZY

  using Distributions,SpecialFunctions,DelimitedFiles,LinearAlgebra,Random,Interpolations
  include("backtr.jl")
  include("create_paramfile.jl")
  include("cova.jl")
  include("krige.jl")
  include("dual.jl")
  include("krige.jl")
  include("nscore.jl")##tableZY
  include("picksupr.jl")
  include("search.jl")
  include("setdual.jl")
  include("setrot.jl")
  include("superblk.jl")
  include("tbmain.jl")
  include("vdc.jl")

    #----------------------------------------------------------------------------------
    # Authors: Xavier Emery and Christian Lantuéjoul, Tom\'as Lagos Gonz\'alez
    # Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
    #         Gaussian random fields via the turning bands method"
    #----------------------------------------------------------------------------------

  function process_input_line(Str)
    splited_line=split(Str, '%');
    j=length(splited_line[1]);
    if j==0 return ""; end
    while splited_line[1][j]==' '
      j=j-1
    end
    return splited_line[1][1:j];
  end
  function smpldist(dist,rng,m,n)
    return reduce(hcat,[ [rand(dist) for 𝔪=1:m] for 𝔫=1:n]);
  end
  function betarnd(α,β,rng,m=1,n=1)
    dist=Beta(α,β);
    return smpldist(dist,rng,m,n);
  end
  function gamrnd(α, θ,rng,m=1,n=1)
    dist=Gamma(α, θ);
    return smpldist(dist,rng,m,n);
  end
  function tbsim(
      α#==2.5=#,simucoord#==Float64[]=#,datacoord=[],ydata=[],nrealiz = 10,return_a_mapping_to_grades=true,nargin=2,nx=200,ny=300,nz=1,
      x0=1.0,y0=1.0,z0=100.0,dx=2.0,dy=2.0,dz=10.0,nd = Int16[1,1,1],limits=[-90,90],
      zmin=0.0,zmax = 10.0,tail=[1.0,5.0],nst=2,nugget=0.1,model=[],
      model1=[1 0.45 100 100 150 0 0 0 1 1000; 2 0.45 100 100 1000000000 0 0 0 1 1000],
      seed = 9784498, radius = [100 100 150] , angles =[0 0 0],
      octant = 1,ndata = 4,nbdecimal=3,ntok = 5000
      )#simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail,
          #model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok)

            # Conditional simulation of a Gaussian random field via the turning bands method
            #-------------------------------------------------------------------------------
            #
            # USE: tbsim(simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail, ...
            #            model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok);
            #
            # INPUT:
            #   simucoord    : coordinates of the locations to simulate (void if a grid simulation is required)
            #   x0,y0,z0     : if grid: minimum grid coordinates along x, y and z directions
            #   nx,ny,nz     :          number of grid nodes along x, y and z directions
            #   dx,dy,dz     :          grid meshes along x, y and z directions
            #   nd           : block discretization along x, y and z directions (1 * 3 vector) ([1 1 1] for point-support simulation)
            #                   The block size is given by the grid meshes dx,dy,dz
            #   datacoord    : data coordinates (n * 3 matrix; void for non-conditional simulations)
            #   ydata        : Gaussian conditioning data (n * 1 vector; void for non-conditional simulations)
            #   limits       : trimming limits (inf and sup) for the Gaussian data (1 * 2 vector)
            #   tableZY      : conversion table between original and Gaussian values (void if no transformation is required)
            #                   The first column of the table contains the original values, the second column their normal scores
            #   zmin,zmax    : minimum and maximum values for the original variable
            #   tail         : parameters lambda and lambda prime for lower-tail and upper-tail extrapolations
            #   model        : covariance model for the Gaussian random field (nst * 7 matrix, where nst is the number of nested structures)
            #                  Each row refers to a nested structure and is codified as: [type, scale factors, angles]
            #                  There are three scale factors (along the rotated y, x and z axes) and three angles
            #                  to define the coordinate rotation (azimuth, dip and plunge), see Deutsch and Journel, 1992, p. 25
            #                  Available types:
            #                    1: spherical
            #                    2: exponential
            #                    3: gamma (parameter b > 0)
            #                    4: stable (parameter b between 0 and 2)
            #                    5: cubic
            #                    6: Gaussian
            #                    7: cardinal sine
            #                    8: Bessel J (parameter b > 0)
            #                    9: Bessel K (parameter b > 0)
            #                   10: generalized Cauchy (parameter b)
            #                   11: exponential sine
            #                   12: linear generalized covariance
            #                   13: power generalized covariance (parameter b > 0)
            #                   14: mixed power generalized covariance (parameter b between 0 and 2)
            #                   15: spline generalized covariance (parameter b = even integer)
            #   cc           : sills or slopes of the nested structures (nst * 1 vector)
            #   b            : third parameters for covariance definition (nst * 1 vector) (used for covariance types 3,4,8,9,10 and 13)
            #   nugget       : nugget effect variance
            #   nlines       : number of lines to use for each nested structure (nst * 1 vector)
            #   nrealiz      : number of realizations to draw
            #   seed         : seed number for generating random values
            #   radius       : maximum search radii along y, x and z (rotated system) for conditioning data
            #                   if radius = inf (infinite), a unique neighborhood is assumed and dual kriging is used
            #   angles       : angles for anisotropic search, according to the GSLIB conventions (Deutsch and Journel, 1992, p. 25)
            #   octant       : divide the neighborhood in octants? 1=yes, 0=no
            #   ndata        : number of conditioning data per octant (if octant=1) or in total (if octant=0)
            #   name         : name of output file
            #   nbdecimal    : number of decimals for the values in output file
            #   header       : create a GSLIB header in the output file? 1=yes, 0=no
            #   ntok         : maximum number of points to keep in memory and simulate simultaneously (optional)
            #                   This number defines how many locations are projected onto the lines at each step of the simulation

            #-------------------------------------------------------------------------------------------------------------------------------
            # This program uses the following subroutines:
            #   backtr.m           : back-transformation from Gaussian to original values
            #   betarnd.m          : simulate beta values
            #   cova.m             : compute covariance values
            #   create_paramfile.m : create default parameter file
            #   DelimitedFiles.readdlm.m          : import ASCII data file
            #   dual.m             : compute dual kriging weights
            #   gamrnd.m           : simulate gamma values
            #   krige.m            : compute simple or ordinary kriging weights
            #   picksupr.m         : build template of super-blocks centered at the block containing the node to simulate
            #   rand.m             : simulate uniform values
            #   randn.m            : simulate normal values
            #   search.m           : search the indices of the data located in the neighborhood
            #   setdual.m          : set up right hand side member for dual kriging
            #   setrot.m           : set up matrix for rotation and reduction of coordinates
            #   superblk.m         : set up super-block strategy
            #   tbmain.m           : main routine for random number generation along the lines
            #   vdc.m              : create equidistributed directions on the sphere (Van der Corput sequence)
            #-----------------------------------------------------------------------------------------------------------

            # Prompt for parameter file if no input is entered
            #-------------------------------------------------

      fid=0;
      if (nargin < 1)

        println("Which parameter file do you want to use?");
        # paramfile = input('','s');
        paramfile="tbsim.par";
      end
      seed=abs(rand(Int64));

      # If a single input is entered, it is the parameter file
      #-------------------------------------------------------

      # if (nargin == 1)
      #   paramfile = simucoord;
      # end


      # Read from parameter file
      #-------------------------

      if (nargin < 2)

        paramfile = "tbsim.par";

        # fid = open(paramfile);
        #
        # if (fid < 0)
        #   close(paramfile);#"all");
        #   println("ERROR - The parameter file does not exist,");
        #   println("        Check for the file and try again");
        #   println(" ");
        #   println("        creating a blank parameter file");
        #   println(" ");
        #   println("Stop - Program terminated.");
        #   create_paramfile;
        #   return;
        # end
        # close(paramfile)
        try
          fid = open(paramfile);
        catch SystemError
          create_paramfile()
          fid = open(paramfile);
        end

        lines_file = readlines(fid);
        close(fid);
        # lines_file=lines_file[2:end];
        option = parse(Int16,process_input_line(lines_file[5]));# type of simulation: 0=gridded locations; 1=scattered locations
        splited_line = process_input_line(lines_file[6]);
        if option ==1
          simucoord = DelimitedFiles.readdlm(splited_line);
        end

        index = parse.(Int64,split(process_input_line(lines_file[7]), ' '));
        if option == 1
          simucoord = simucoord[:,index];
        else
          simucoord = Float64[];
        end

        # # The parameter file does exist
        # fgets(fid); fgets(fid); fgets(fid); fgets(fid);
        #
        # tline = fgets(fid);
        # option = str2num(tline);# type of simulation: 0=gridded locations; 1=scattered locations
        #
        # tline = fgets(fid);
        # i = (find(tline == " "));
        # if ~isempty(i), tline = tline(1:i-1); end
        # fid2 = fopen(tline);
        # if (option == 1) # scattered locations
        #   if (fid2 < 0)
        #     fclose("all");
        #     error(["Cannot import file ",tline]);
        #   else
        #     simucoord = DelimitedFiles.readdlm(tline);
        #   end
        # end

        # tline = fgets(fid);
        # index = str2num(tline);
        # if (option == 1), simucoord = simucoord(:,index); else, simucoord = []; end

        if option == 0
          splited_line = split(process_input_line(lines_file[8]), ' ');
          x0=parse(Float16,splited_line[1]);
          y0=parse(Float16,splited_line[2]);
          z0=parse(Float16,splited_line[3]);
          splited_line = split(process_input_line(lines_file[9]), ' ');
          nx=parse(Int16,splited_line[1]);
          ny=parse(Int16,splited_line[2]);
          nz=parse(Int16,splited_line[3]);
        end
        splited_line = split(process_input_line(lines_file[10]), ' ');
        dx = parse(Float16,splited_line[1]);
        dy = parse(Float16,splited_line[2]);
        dz = parse(Float16,splited_line[3]);

        splited_line = split(process_input_line(lines_file[11]), ' ');
        nd = [parse(Int16,splited_line[1]),
              parse(Int16,splited_line[2]),
              parse(Int16,splited_line[3])];

        # tline = fgets(fid);
        # tline = str2num(tline);
        # if (option == 0), x0 = tline(1); y0 = tline(2); z0 = tline(3); end

        # tline = fgets(fid);
        # tline = str2num(tline);
        # if (option == 0), nx = tline(1); ny = tline(2); nz = tline(3); end

        # tline = fgets(fid);
        # tline = str2num(tline);
        # dx = tline(1); dy = tline(2); dz = tline(3);

        # tline = fgets(fid);
        # nd = str2num(tline);

        splited_line = process_input_line(lines_file[12]);
        ydata = DelimitedFiles.readdlm(splited_line);

        # tline = fgets(fid);
        # i = (find(tline == " "));
        # if ~isempty(i), tline = tline(1:i-1); end
        # fid3 = fopen(tline);
        # if (fid3 > -1), ydata = DelimitedFiles.readdlm(tline); end

        splited_line = split(process_input_line(lines_file[13]), ' ');
        if length(splited_line)>0
          index=[parse(Int16,splited_line[1]),
                parse(Int16,splited_line[2]),
                parse(Int16,splited_line[3])];
          datacoord = ydata[:,index];
        else
          datacoord = Float64[];
        end

        # tline = fgets(fid);
        # if (fid3 > -1)
        #   index = str2num(tline);
        #   datacoord = ydata(:,index);
        # else
        #   datacoord = [];
        # end

        splited_line = process_input_line(lines_file[14]);
        if length(splited_line)>0
          ydata = ydata[:,parse(Int16,splited_line)];
        else
          ydata=Float64[];
        end

        # tline = fgets(fid);
        # if (fid3 > -1)
        #   index = str2num(tline);
        #   ydata = ydata(:,index);
        # else
        #   ydata = [];
        # end

        splited_line = split(process_input_line(lines_file[15]), "  ");
        limits = [parse(Int16,splited_line[1]),parse(Int16,splited_line[2])];

        # tline = fgets(fid);
        # limits = str2num(tline);

        # splited_line = process_input_line(lines_file[16]);
        # if length(splited_line)>0
        #   tableZY = readdlm(splited_line);
        # end

        # tline = fgets(fid);
        # i = (find(tline == " "));
        # if ~isempty(i), tline = tline(1:i-1); end
        # fid4 = fopen(tline);
        # if (fid4 > -1), tableZY = DelimitedFiles.readdlm(tline); else, tableZY = []; end

        splited_line = split(process_input_line(lines_file[17]), ' ');
        zmin,zmax = parse(Float16,splited_line[1]),parse(Float16,splited_line[2]);

        splited_line = split(process_input_line(lines_file[18]), ' ');
        tail = [parse(Float16,splited_line[1]),parse(Float16,splited_line[2])];

        splited_line = split(process_input_line(lines_file[19]), ' ');
        nst,nugget = parse(Int16,splited_line[1]),parse(Float16,splited_line[2]);

        # tline = fgets(fid);
        # tline = str2num(tline);
        # zmin = tline(1); zmax = tline(2);

        # tline = fgets(fid);
        # tail = str2num(tline);

        # tline = fgets(fid);
        # tline = str2num(tline);
        # nst = tline(1); nugget = tline(2);

        model=[];
        cc=Array{Float16}(undef,nst,1);
        b=Array{Int16}(undef,nst,1);
        nlines=Array{Int16}(undef,nst,1);
        for i = 1:nst
          splited_line = split(process_input_line(lines_file[19+i]), ' ');
          dd=[ j==1 ? parse(Float32,splited_line[j]) : parse(Int32,splited_line[j]) for j in [1 3 4 5 6 7 8] ];
          model=vcat(model,dd);
          cc[i,1] = parse(Float16,splited_line[2]);
          b[i,1] = parse(Int16,splited_line[9]);
          nlines[i,1] = parse(Int16,splited_line[10]);

          # tline = fgets(fid);
          # tline = str2num(tline);
          # model(i,:) = tline([1 3 4 5 6 7 8]);
          # cc(i,1) = tline(2);
          # b(i,1) = tline(9);
          # nlines(i,1) = tline(10);
        end

        nrealiz = parse(Int16,process_input_line(lines_file[20+nst]));
        seed = parse(Int32,process_input_line(lines_file[21+nst]));#abs(rand(Int32));#
        radius = [ parse(Int16,aux) for aux in split(process_input_line(lines_file[22+nst]), ' ')];
        angles = [ parse(Float16,aux) for aux in split(process_input_line(lines_file[23+nst]), ' ')];
        octant = parse(Int16,process_input_line(lines_file[24+nst]));
        ndata = parse(Int16,process_input_line(lines_file[25+nst]));
        name = process_input_line(lines_file[26+nst])
        nbdecimal = parse(Int16,process_input_line(lines_file[27+nst]));
        header = parse(Int16,process_input_line(lines_file[28+nst]));
        ntok = parse(Int32,process_input_line(lines_file[29+nst]));

        # tline = fgets(fid);
        # nrealiz = str2num(tline);
        #
        # tline = fgets(fid);
        # seed = str2num(tline);
        #
        # tline = fgets(fid);
        # radius = str2num(tline);
        #
        # tline = fgets(fid);
        # angles = str2num(tline);
        #
        # tline = fgets(fid);
        # octant = str2num(tline);
        #
        # tline = fgets(fid);
        # ndata = str2num(tline);
        #
        # name = fgets(fid);
        # i = (find(name == " "));
        # if ~isempty(i), name = name(1:i-1); end
        #
        # tline = fgets(fid);
        # nbdecimal = str2num(tline);
        #
        # tline = fgets(fid);
        # header = str2num(tline);
        #
        # tline = fgets(fid);
        # ntok = str2num(tline);
        #
        # fclose("all");

      end


      # Define default values
      #----------------------
      if isempty(model)
        nst = size(model1,1); # number of nested structures
        cc=Array{Float16}(undef,nst,1);
        b=Array{Int16}(undef,nst,1);
        nlines=Array{Int16}(undef,nst,1);
        for i = 1:nst
          if i==1
            model=[ model1[i,j] for j in [1 3 4 5 6 7 8] ];
          else
            dd=[ model1[i,j] for j in [1 3 4 5 6 7 8] ];
            model=vcat(model,dd);
          end
          cc[i,1] = Float16(model1[i,2]);
          b[i,1] = Int16(model1[i,9]);
          nlines[i,1] = Int16(model1[i,10]);
        end
      end
      # warning("off","all");
      max_intervalnumber = 1e4;
      block = [dx dy dz];  # block dimensions
      ng = prod(nd);       # number of points discretizing the block
      sigma_nugget = sqrt(nugget);
      # cc = cc(:);
      # b = b(:);
      # nlines = nlines(:);
      if ng == 1
        block = [0 0 0];
      end


      # Define the grid that discretizes the block
      #-------------------------------------------

      t2 = Array{Float64}(undef,ng,3);
      for i = 1:3
        nl = prod(nd[1:i-1]);
        nr = prod(nd[i+1:3]);
        t1  = reshape([0.5*(1/nd[i]-1):1/nd[i]:0.5*(1-1/nd[i])][1],(nd[i],1));
        t2[:,i]=kron(ones(nl,1),kron(t1,ones(nr,1)));
      end
      grid = t2.*(ones(ng,1)*block);

      # println("grid=",grid);
      # Check input
      #------------

      # if isempty(ntok), ntok = 5000; end

      if length(nlines) != nst
        nlines = nlines[1]*ones(nst,1);
      end

      index = Int16[];
      for i = 1:nst
        if cc[i] > eps() # valid nested structure
          push!(index,i);
        end
      end
      cc = cc[index];
      model = model[index,:];
      b = b[index];
      nlines = nlines[index];
      nst = size(model,1);
      order = -ones(nst,1);

      if maximum(model[:,1]) < 11.5 # stationary model
        flag = 0;
        variance = round(sum(cc)+nugget,digits=3);
        if abs(variance - 1) > eps()
          println(" ");
          println("WARNING - The variogram sill is not equal to 1");
          println(" ");
        end
      else # unbounded variogram model
        flag = 1;
      end

      if length(radius) != 3
        radius = radius[1]*ones(1,3);
      end

      # if isempty(tableZY) && length(ydata)>0
      #   # tableZY=hcat(sort(exp.(ydata/α)),sort(randn(length(ydata))));
      #   tableZY=readdlm("nscore.trn");
      #   if !isempty(tableZY)
      #     # tableZY = tableZY[:,1:2];
      #     length_table = size(tableZY,1);
      #     tableZY[:,2] = tableZY[:,2] + 1e-10*ones(length_table); # avoid tied values in the conversion table
      #   end
      # end

      for i = 1:nst

        if (model[i,1] == 3) # Gamma model
          if (b[i] <= 0)
            error("The parameter of the gamma model must be positive");
          end

        elseif (model[i,1] == 4) # Stable model
          if (b[i] > 2)
            error("The parameter of the stable model must be less than 2");
          elseif (b[i] <= 0)
            error("The parameter of the stable model must be positive");
          elseif (b[i] == 2) # Gaussian model
            model[i,1] = 6;
          elseif (b[i] == 1) # Exponential model
            model[i,1] = 2;
          elseif (b[i] > 1) # Simulation via spectral method
            model[i,1] = 4.5;
          end

        elseif (model[i,1] == 8) # Bessel-J model
          if (b[i] < 0.5)
            error("The parameter of the Bessel-J model must be greater than 0.5");
          elseif (b[i] == 0.5) # Cardinal sine model
            model[i,1] = 7;
          end

        elseif (model[i,1] == 9) # Bessel-K model
          if (b[i] <= 0)
            error("The parameter of the Bessel-K model must be positive");
          elseif (b[i] == 0.5) # Exponential model
            model[i,1] = 2;
          elseif (b[i] > 0.5) # Simulation via spectral method
            model[i,1] = 9.5;
          end

        elseif (model[i,1] == 10) # Generalized Cauchy model
          if (b[i] <= 0)
            error("The parameter of the generalized Cauchy model must be positive");
          end

        elseif (model[i,1] == 12) # Linear model
          order[i] = 0;

        elseif (model[i,1] == 13) # Power model
          if (b[i] <= 0)
            error("The exponent of the power model must be positive");
          end
          if (b[i]/2 == floor(b[i]/2))
            model[i,1] = 13.1;
          elseif (b[i] > 1)
            model[i,1] = 13.5;
          end
          order[i] = ceil(b[i]/2)-1;

        elseif (model[i,1] == 14) # Mixed power model
          if (b[i] <= 0)
            error("The exponent of the mixed power model must be positive");
          end
          if (b[i] > 2)
            error("The exponent of the mixed power model must be less than 2");
          end
          order[i] = 0;

        elseif (model[i,1] == 15) # Spline model
          if (b[i] < 0)
            error("The exponent of the spline model must be positive");
          end
          if (b[i]/2 != floor(b[i]/2))
            error("The exponent of the spline model must an even integer");
          end
          order[i] = b[i]/2;
        end

      end

      maxorder = maximum(order);


      # Assign conventional values to the grid parameters (x0,y0,z0,nx,ny,nz,dx,dy,dz)
      # if the simulation is to be performed on scattered locations
      #-------------------------------------------------------------------------------

      if !isempty(simucoord)
        if size(simucoord,1) > 1
          minsimucoord = minimum(simucoord,dims=1);
          maxsimucoord = maximum(simucoord,dims=1);
        else
          minsimucoord = simucoord;
          maxsimucoord = simucoord;
        end
        x0 = minsimucoord[1];
        y0 = minsimucoord[2];
        z0 = minsimucoord[3];
        if maxsimucoord[1]>minsimucoord[1]+eps()
          nx = 100;
          dx = (maxsimucoord[1]-minsimucoord[1])/(nx-1);
        else
          nx = 1;
          dx = 1;
        end
        if maxsimucoord[2]>minsimucoord[2]+eps()
          ny = 100;
          dy = (maxsimucoord[2]-minsimucoord[2])/(ny-1);
        else
          ny = 1;
          dy = 1;
        end
        if maxsimucoord[3]>minsimucoord[3]+eps()
          nz = 100;
          dz = (maxsimucoord[3]-minsimucoord[3])/(nz-1);
        else
          nz = 1;
          dz = 1;
        end
      end


      # Remove data whose values are not in the trimming limits interval
      #-----------------------------------------------------------------

      m0 = size(datacoord,1);

      if m0 > 0
        I=(LinearIndices(ydata))[findall( xy->(xy > limits[1] && xy<limits[2]), ydata)];
        # I = find( (ydata > limits(1)) & (ydata < limits(2)) );
        datacoord = datacoord[I,:];
        ydata = ydata[I];
      # end
      # println("mean(ydata)=",mean(ydata))
      # println("std(ydata)=",std(ydata))

      # Remove data located too far from the locations to simulate
      #-----------------------------------------------------------

      # m0 = size(datacoord,1);
      #
      # if m0 > 0

        flag = 2*flag;

        # Reduced-rotated data coordinates
        search_rotationmatrix = setrot(vcat([1],radius',angles')',1); # rotation-reduction matrix for data search
        tmp = datacoord*search_rotationmatrix;

        # Extremal points to simulate
        x = [x0 y0 z0;
             x0 y0 z0+(nz-1)*dz;
             x0 y0+(ny-1)*dy z0;
             x0 y0+(ny-1)*dy z0+(nz-1)*dz;
             x0+(nx-1)*dx y0 z0;
             x0+(nx-1)*dx y0 z0+(nz-1)*dz;
             x0+(nx-1)*dx y0+(ny-1)*dy z0;
             x0+(nx-1)*dx y0+(ny-1)*dy z0+(nz-1)*dz];
        x = x*search_rotationmatrix;
        minx = minimum(x[:,1]);
        miny = minimum(x[:,2]);
        minz = minimum(x[:,3]);
        maxx = maximum(x[:,1]);
        maxy = maximum(x[:,2]);
        maxz = maximum(x[:,3]);

        # Identify and remove the data located beyond the search radii
        I = findall( xy->( (xy[1] < minx-1) || (xy[1] > maxx+1) ||
              (xy[2] < miny-1) || (xy[2] > maxy+1) ||
              (xy[3] < minz-1) || (xy[3] > maxz+1) ), [tmp[xx,:] for xx=1:length(tmp[:,1])]);
            # find( (tmp(:,1) < minx-1) | (tmp(:,1) > maxx+1) |
            #       (tmp(:,2) < miny-1) | (tmp(:,2) > maxy+1) |
            #       (tmp(:,3) < minz-1) | (tmp(:,3) > maxz+1) );
        datacoord = datacoord[setdiff(1:length(datacoord[:,1]),I),:] ;
        ydata = ydata[setdiff(1:length(ydata[:,1]),I),:] ;
        m0 = size(datacoord,1);

        # println("mean(ydata)=",mean(ydata))
        # println("std(ydata)=",std(ydata))
      end


      # Create seed numbers
      #--------------------

      # rand("state",seed);
      # randn("state",seed);
      rng = Random.MersenneTwister(seed);
      seed_vdc = abs.(Int.(rand(rng,Int64,nst)));
      seed_line = abs.(Int.(rand(rng,Int64,nst*nrealiz,maximum(nlines))));

      # Extremal coordinates to simulate (including the data locations) in the original referential
      #--------------------------------------------------------------------------------------------

      if m0 > 0

        minx = minimum(vcat(datacoord[:,1],[x0-block[1]/2]));
        miny = minimum(vcat(datacoord[:,2],[y0-block[2]/2]));
        minz = minimum(vcat(datacoord[:,3],[z0-block[3]/2]));
        maxx = maximum(vcat(datacoord[:,1],[x0+(nx-1)*dx+block[1]/2]));
        maxy = maximum(vcat(datacoord[:,2],[y0+(ny-1)*dy+block[2]/2]));
        maxz = maximum(vcat(datacoord[:,3],[z0+(nz-1)*dz+block[3]/2]));

      else

        radius = [0 0 0];
        minx = x0-block[1]/2;
        miny = y0-block[2]/2;
        minz = z0-block[3]/2;
        maxx = x0+(nx-1)*dx+block[1]/2;
        maxy = y0+(ny-1)*dy+block[2]/2;
        maxz = z0+(nz-1)*dz+block[3]/2;

      end

      extreme_coord = [minx miny minz; minx miny maxz; minx maxy minz; minx maxy maxz;
                       maxx miny minz; maxx miny maxz; maxx maxy minz; maxx maxy maxz];


      #--------------------------------------------------------------------------------------------

      # PREPARE THE LINES
      #------------------

      # Initialization
      all_lines = Array{Float64}(undef,nrealiz*maximum(nlines),3,nst);
      all_offset = Array{Float64}(undef,nrealiz*maximum(nlines),nst);
      all_r = Array{Float64}(undef,maximum(nlines),nrealiz,nst);
      all_phi = Array{Float64}(undef,maximum(nlines),nrealiz,nst);
      all_theta = Array{Float64}(undef,maximum(nlines),nrealiz,nst);
      valid_lines = Array{Float64}(undef,nrealiz*maximum(nlines),nst);
      total_nugget = nugget*ones(1,nrealiz);
      sigma=zeros(nst);

      model_rotationmatrix=Array{Float32}(undef,3,3,nst);

      for i = 1:nst

        # Line generation (Van der Corput sequence)
        #------------------------------------------
        V = vdc(nlines[i],nrealiz,Int(seed_vdc[i]));
        R = setrot(model,i);
        # println("mean(R)",mean(R,1))
        model_rotationmatrix[:,:,i] = R;
        lines = V*R';
        all_lines[1:nrealiz*nlines[i],:,i] = lines;
        # println("mean(lines)",mean(lines,1))


        # Spherical covariance model?
        #----------------------------
        if (model[i,1] == 1)
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(3*cc[i]./nlines[i]);
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)'-minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,(nlines[i],nrealiz));
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Exponential covariance model?
        #------------------------------
      elseif (model[i,1] == 2)
          v = randn(rng,6,nlines[i]*nrealiz);
          w = 3*rand(rng,nlines[i]*nrealiz);
          v[5,:] = v[5,:].*(w.>1);
          v[6,:] = v[6,:].*(w.>1);
          G = 0.5*sum(v.^2)' * ones(1,3);
          lines = lines./G;
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(3*cc[i]./nlines[i]);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)' - minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,(nlines[i],nrealiz));
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Gamma covariance model?
        #------------------------
      elseif (model[i,1] == 3)
          v = randn(rng,6,nlines[i]*nrealiz);
          w = 3*rand(rng,nlines[i]*nrealiz);
          v[5,:] = v[5,:].*(w.>1);
          v[6,:] = v[6,:].*(w.>1);
          G = 0.5*(sum(v.^2)'./gamrnd(b[i],1,rng,nlines[i]*nrealiz,1)) * ones(1,3);
          lines = lines./G;
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(3*cc[i]./nlines[i]);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)' - minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,(nlines[i],nrealiz));
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Stable covariance model with parameter <1?
        #-------------------------------------------
        elseif (model[i,1] == 4)
          v = randn(rng,6,nlines[i]*nrealiz);
          w = 3*rand(rng,1,nlines[i]*nrealiz);
          v[5,:] = v[5,:].*(w>1);
          v[6,:] = v[6,:].*(w>1);
          e = -log(rand(rng,nlines[i]*nrealiz,1));
          f = pi*rand(rng,nlines[i]*nrealiz,1)-pi/2;
          stable = abs( sin(b[i]*(f-pi/2))./(cos(f).^(1 ./ b[i])).*(cos(f-b[i]*(f-pi/2))./e).^((1-b[i])./b[i]) );
          G = 0.5*(sum(v.^2)'./stable) * ones(1,3);
          lines = lines./G;
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(3*cc[i]./nlines[i]);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)'-minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,(nlines[i],nrealiz));
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Stable model with parameter >1?
        #--------------------------------
        elseif (model[i,1] == 4.5)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2;
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(6*t);
          e = -log.(rand(rng,nlines[i]*nrealiz,1));
          f = pi*rand(rng,nlines[i]*nrealiz,1) .- pi/2;
          G = abs.( sin.(b[i]/2*(f .- pi/2))./(cos.(f).^(2 ./ b[i])).*(cos.(f .- b[i]/2*(f .- pi/2)) ./ e).^((1-b[i]/2) ./ b[i]*2) )*ones(1,3);
          lines = lines.*sqrt.(G/3);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Cubic covariance model?
        #------------------------
        elseif (model[i,1] == 5)
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(210*cc[i]./nlines[i]);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)'-minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,(nlines[i],nrealiz));
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Gaussian covariance model?
        #---------------------------
        elseif (model[i,1] == 6)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2;
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(2*t);
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Cardinal sine covariance model?
        #--------------------------------
        elseif (model[i,1] == 7)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = 2*rand(rng,nlines[i],nrealiz) .- 1;
          all_r[1:nlines[i],1:nrealiz,i] = sign.(t);
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Bessel-J covariance model?
        #---------------------------
        elseif (model[i,1] == 8)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = betarnd(1.5,b[i]-0.5,rng,nlines[i],nrealiz);
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(t);
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Bessel-K covariance model (b<0.5)?
        #-----------------------------------
        elseif (model[i,1] == 9)
          v = randn(rng,6,nlines[i]*nrealiz);
          w = 3*rand(rng,1,nlines[i]*nrealiz);
          v[5,:] = v[5,:].*(w>1);
          v[6,:] = v[6,:].*(w>1);
          G = 0.5*(sum(v.^2)'.*sqrt(betarnd(b[i],0.5-b[i],rng,nlines[i]*nrealiz,1))) * ones(1,3);
          lines = lines./G;
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt(3*cc[i]./nlines[i]);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);

          # Identify the lines for which the scale factor is too small
          # The simulation along these lines will be replaced by a nugget effect
          interval_number = maximum(x,dims=1)' - minimum(x,dims=1)';
          nottoosmall = (interval_number .< max_intervalnumber);
          valid_lines[1:nrealiz*nlines[i],i] = nottoosmall;
          invalid_lines = reshape(1 .- nottoosmall,nlines[i],nrealiz);
          total_nugget = total_nugget .+ sum(invalid_lines)*cc[i]./nlines[i];


        # Bessel-K covariance model (b>0.5)?
        #-----------------------------------
        elseif (model[i,1] == 9.5)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2;
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(6*t);
          G = gamrnd(b[i],1,rng,nlines[i]*nrealiz,1)*ones(1,3);
          lines = lines ./ sqrt.(G*12);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Generalized Cauchy model?
        #--------------------------
        elseif (model[i,1] == 10)
          sigma[i] = sqrt(2*cc[i]./nlines[i]);
          t = randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2;
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(6*t);
          G = gamrnd(b[i],1,rng,nlines[i]*nrealiz,1)*ones(1,3);
          lines = lines .* sqrt.(G/3);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Exponential sine model?
        #------------------------
        elseif (model[i,1] == 11)
          t2=Array{Float64}(undef,nlines[i],nrealiz);
          sigma[i] = sqrt.(2*cc[i]./nlines[i]);
          for k = 1:nrealiz
            for j = 1:nlines[i]
              iflag = -1;
              while (iflag < 0)
                t2[j,k] = rand(Gamma(0.5,1));
                u = rand(rng);
                c = 1;
                puis = 1;
                n = 0;
                iflag = 0;
                while (iflag == 0)
                  n = n+1;
                  puis = puis * pi*pi/4 / (2*n*(2*n-1));
                  if (n/2 == round(n/2))
                    c = c+puis*exp(-4*t2[j,k]*n*(n+1));
                    if (u > c) iflag = -1; end
                  else
                    c = c-puis*exp(-4*t2[j,k]*n*(n+1));
                    if (u < c) iflag = 1; end
                  end
                end
              end
            end
          end
          t1 = randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2 + randn(rng,nlines[i],nrealiz).^2;
          all_r[1:nlines[i],1:nrealiz,i] = sqrt.(t1 ./ 2 ./ t2);
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);


        # Linear model?
        #--------------
        elseif (model[i,1] == 12)
          x = extreme_coord*lines';
          sigma[i] = 2*sqrt.(cc[i]./nlines[i]);

          # Modify the slope and scale factor if the latter is too small
          delta = maximum(x,dims=1) - minimum(x,dims=1);
          factor = maximum(delta)/max_intervalnumber;
          if (factor > 1)
            lines = lines./factor;
            sigma[i] = sigma[i]*sqrt(factor);
            x = extreme_coord*lines';
          end

          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);


        # Power model with exponent <=1?
        #-------------------------------
        elseif (model[i,1] == 13)

          # Set the scale factor equal to the diameter of the domain to simulate
          x = extreme_coord*lines';
          delta = maximum(x,dims=1)-minimum(x,dims=1);
          factor = maximum(delta);
          model[i,2:4] = factor*model[i,2:4];
          cc[i] = (factor.^b[i])*cc[i];
          R = setrot(model,i);
          model_rotationmatrix[:,:,i] = R;
          lines = V*R';

          # The power covariance is simulated as a mixture of triangular models
          sigma[i] = sqrt(2*gamma(b[i]/2+1.5)/gamma(b[i]/2+0.5)*cc[i]./nlines[i]);
          u = rand(rng,nlines[i]*nrealiz,1);
          G = max.(0.001,min.(1,(u/(1-b[i])).^(1/b[i]))) * ones(1,3);
          lines = lines./G;
          x = extreme_coord*lines';

          # Modify the slope and scale factor if the latter is too small
          delta = maximum(x,dims=1) - minimum(x,dims=1);
          factor = maximum(delta)/max_intervalnumber;
          if (factor > 1)
            model[i,2:4] = factor*model[i,2:4];
            cc[i] = (factor.^b[i])*cc[i];
            R = setrot(model,i);
            model_rotationmatrix[:,:,i] = R;
            lines = V*R';
            lines = lines./G;
            sigma[i] = sqrt(2*gamma(b[i]/2+1.5)/gamma(b[i]/2+0.5)*cc[i]./nlines[i]);
            x = extreme_coord*lines';
          end

          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_offset[1:nrealiz*nlines[i],i] = minimum(x,dims=1)' - rand(rng,nrealiz*nlines[i],1);


        # Pure drift model?
        #------------------
        elseif (model[i,1] == 13.1)
          sigma[i] = sqrt(cc[i]./nlines[i]);
          all_r[1:nlines[i],1:nrealiz,i] = exp(0.5*lgamma(2*order[i]+4)-lgamma(order[i]+2))*randn(rng,nlines[i],nrealiz);
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_phi[1,1,i] = order[i]+1;


        # Power model with exponent >1?
        #------------------------------
        # elseif (model[i,1] == 13.5)
        #   sigma[i] = sqrt(cc[i]./nlines[i]*gamma(b[i]+2)/gamma(b[i]+2-2*order[i]));
        #   t = (randn(nlines[i],nrealiz)./randn(nlines[i],nrealiz)).^2;
        #   all_r[1:nlines[i],1:nrealiz,i] = 2*pi*t;
        #   all_lines[1:nrealiz*nlines[i],:,i] = lines;
        #   all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(nlines[i],nrealiz);
        #   all_theta[1:nlines[i],1:nrealiz,i] = sqrt.( (-pi^(0.5-b[i])*gamma(b[i]/2+1.5-order[i])/
        #             2^(2*order[i]-3)/gamma(order[i]-b[i]/2)) .* ( (1+t)./(t.^(b[i]+0.5)))) ;


        # Mixed power model?
        #-------------------
        elseif (model[i,1] == 14)
          sigma[i] = sqrt(b[i]*cc[i]./nlines[i]);
          t = (randn(rng,nlines[i],nrealiz) ./ randn(rng,nlines[i],nrealiz)).^2;
          all_r[1:nlines[i],1:nrealiz,i] = 2*pi*t;
          all_lines[1:nrealiz*nlines[i],:,i] = lines;
          all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(rng,nlines[i],nrealiz);
          pow = rand(rng,nlines[i],nrealiz)*b[i];
          all_theta[1:nlines[i],1:nrealiz,i] = sqrt.( -8 * pi .^ (0.5 .- pow) .* gamma.(pow/2 .+ 1.5) ./
                    gamma.(-pow/2).*(1 .+ t)./exp.((pow .+ 0.5).*log.(t)) );


        # Spline model?
        #--------------
        # elseif (model[i,1] == 15)
        #   sigma[i] = sqrt(cc[i]./nlines[i]);
        #   t = (randn(nlines[i],nrealiz)./randn(nlines[i],nrealiz)).^2;
        #   all_r[1:nlines[i],1:nrealiz,i] = 2*pi*t;
        #   all_lines[1:nrealiz*nlines[i],:,i] = lines;
        #   all_phi[1:nlines[i],1:nrealiz,i] = 2*pi*rand(nlines[i],nrealiz);
        #   all_theta[1:nlines[i],1:nrealiz,i] = sqrt( pi.^(1-b[i]).*gamma(b[i]+2)./
        #             2.^(b[i]-1).*(1+t)./t.^(b[i]+0.5) );

        end

      end

      total_nugget = sqrt.(total_nugget);
      max_nugget = maximum(total_nugget);


      #--------------------------------------------------------------------------------------------

      # NON CONDITIONAL SIMULATION AT THE DATA LOCATIONS
      #-------------------------------------------------

      simudata = zeros(m0,nrealiz);
      weights = 0;
      outputformat = [];
      if (m0 > 0)

        println(" ");
        println("Simulation at data locations");


        # How many data locations can be simulated simultaneously?
        #---------------------------------------------------------

        m2 = max(1,min(m0,ntok));
        sequence = [ λ for λ = 0:m2:(m0-0.5)];
        push!(sequence,m0);
        seed_nugget_data = abs(Int(rand(rng,Int64)));


        # Loop over the sequences of data points
        #---------------------------------------

        for n = 1:length(sequence)-1

          index = Int32[ λ for λ=(sequence[n]+1):(sequence[n+1])];
          simudata[index,:] = tbmain(datacoord[index,:],model,sigma,nlines,nrealiz,seed_line,
                                     all_lines,all_offset,all_r,all_phi,all_theta,valid_lines);

        end


        # Add nugget effect
        #------------------

        if max_nugget > eps()
          rng=Random.MersenneTwister(seed_nugget_data);
          # randn('state',seed_nugget_data);
          simudata = simudata + (ones(m0,1)*total_nugget).*randn(rng,m0,nrealiz);
        end


        # Prepare conditioning kriging
        #-----------------------------

        if isinf(radius[1]) # kriging in a unique neighborhood

          residuals = ydata*ones(1,nrealiz)-simudata;
          weights = dual(datacoord,residuals,model,cc,b,nugget+1e-7,model_rotationmatrix,maxorder);

        else # prepare super-block search strategy

          I,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup = superblk(datacoord,nx,ny,nz,x0,y0,z0,dx,dy,dz);

          # sort data according to ascending super-block number
          datacoord = datacoord[I,:];
          ydata = ydata[I,:];
          simudata = simudata[I,:];
          residuals = ydata*ones(1,nrealiz)-simudata;

          # build template of super-blocks centered at the block containing the node to simulate
          ixsbtosr,iysbtosr,izsbtosr = picksupr(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,search_rotationmatrix);

        end

      end


      #--------------------------------------------------------------------------------------------

      # PREPARE OUTPUT FILE
      #--------------------

      # Open output file and prepare output format
      #-------------------------------------------

      # fid = fopen(name,"w");
      # outputformat = [];


      # Prepare format for the values in the output file

      # if isempty(tableZY)
      #   nbdigit = nbdecimal[1]+3;
      # else
      #   nbdigit = maximum(abs.(tableZY[:,1]));
      #   nbdigit = floor(log10(nbdigit))+nbdecimal[1]+3;
      # end
      # for i = 1:nrealiz
      #   outputformat = [outputformat," #",int2str(nbdigit),".",int2str(nbdecimal(1)),"f"];
      # end


      # Create a header?

      # if header > 0
      #
      #   fprintf(fid,"#1s\n","Simulated values");
      #   fprintf(fid,"#1s\n",int2str(nrealiz));
      #
      #   for i = 1:nrealiz
      #     fprintf(fid,"#1s\n",["Realization no.",int2str(i)]);
      #   end
      #
      # end

      progress = [0];
      COORDS=[];
      RR=[];
      #--------------------------------------------------------------------------------------------

      # CONDITIONAL SIMULATION AT SCATTERED LOCATIONS
      #----------------------------------------------
      simu=[];
      if !isempty(simucoord)

        println(" ");
        println("Simulation at scattered locations");
        println(" ");


        # How many locations can be simulated simultaneously?
        #----------------------------------------------------

        m1 = size(simucoord,1);
        ntok = ceil(ntok/ng);
        m2 = max(1,min(m1,ntok));
        # sequence  = [[0:m2:m1-0.5] m1];
        sequence = [ λ for λ = 0:m2:(m1-0.5)];
        push!(sequence,m1);
        lengthsequence = length(sequence)-1;
        seed_nugget = abs.(rand(rng,Int64,lengthsequence));


        # Loop over the sequences of points to simulate
        #----------------------------------------------

        for n = 1:lengthsequence

          # Coordinates of the points to simulate
          #--------------------------------------

          index = Int32[ λ for λ=(sequence[n]+1):(sequence[n+1])];
          m1 = length(index);
          coord = kron(simucoord[index,:],ones(ng))-kron(ones(m1,1),grid);
          m2 = m1*ng;


          # Non-conditional simulation
          #---------------------------

          simu = tbmain(coord,model,sigma,nlines,nrealiz,seed_line,all_lines,all_offset,
                        all_r,all_phi,all_theta,valid_lines);


          # Add nugget effect
          #------------------

          if max_nugget > eps()
            # randn('state',seed_nugget(n));
            rng = Random.MersenneTwister(seed_nugget[n]);
            simu = simu + (ones(m2,1)*total_nugget).*randn(rng,m2,nrealiz);
          end


          # Conditioning
          #-------------

          if isinf(radius[1]) # dual kriging

            for i = 1:m2

              # Substitution of residuals
              k0 = setdual(model,coord[i,:],cc,b,datacoord,model_rotationmatrix,maxorder);
              simu[i,:] = simu[i,:] + k0*weights;

            end


          elseif radius[1] > eps() # kriging in a moving neighborhood

            for i = 1:m1

              # Coordinates of the centre of the ith block to simulate
              block_centre = [mean(coord[(i-1)*ng+1:i*ng,1]) mean(coord[(i-1)*ng+1:i*ng,2]) mean(coord[(i-1)*ng+1:i*ng,3])];

              # Search for neighboring data
              datacoord_i,residuals_i = Search(datacoord,residuals,block_centre,search_rotationmatrix,octant,ndata,nxsup,
                               nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb);

              # Substitution of residuals
              if !isempty(datacoord_i) # samples are found in the kriging neighborhood
                weights = krige(datacoord_i,coord[(i-1)*ng+1:i*ng,:],model,cc,b,nugget+1e-7,model_rotationmatrix,maxorder);
                simu[(i-1)*ng+1:i*ng,:] = simu[(i-1)*ng+1:i*ng,:] + weights'*residuals_i;
              elseif (flag > 1) # intrinsic model, no conditioning data found
                println(" ");
                println("WARNING - No data found in moving neighborhood. Corresponding results may be senseless");
                println(" ");
                flag = flag - 1;
              end

            end

          end


          # Back transform to original scale
          #---------------------------------
          # open("share/simu.txt", "w") do ff
          #   writedlm(ff,simu);
          # end
          # open("share/tableZY.txt", "w") do ff
          #   writedlm(ff,tableZY);
          # end
          # open("share/zmin-zmax.txt", "w") do ff
          #   writedlm(ff,[zmin zmax]);
          # end
          # open("share/tail.txt", "w") do ff
          #   writedlm(ff,[tail]);
          # end
          # simu=readdlm("share/simu.txt");
          # tableZY=readdlm("share/tableZY.txt");
          # zmin,zmax=readdlm("share/zmin-zmax.txt");
          # tail=readdlm("share/tail.txt");
          simu = backtr(simu,tableZY,zmin,zmax,tail);
          # println("size(simu)=",size(simu))

          # Average the values over the blocks
          #-----------------------------------

          if ng > 1
            simu = reshape(simu,(ng,m1*nrealiz));
            simu = mean(simu);
            simu = reshape(simu,(m1,nrealiz));
          end
          # println("size(simu)=",size(simu))


          # Report on progress from time to time
          #-------------------------------------

          push!(progress,10*floor((10*n)/length(sequence)));
          if (progress[n+1] > progress[n])
            # disp(['  ',num2str(progress(n+1)),'# completed']);
            sleep(0.001);
          end


          # Write in output file
          #---------------------

          # fprintf(fid,[outputformat,"\n"],simu');
          # println("size(outputformat)=",size(outputformat))
          if n==1
              outputformat=simu
          else
              outputformat=vcat(outputformat,simu);
          end

        end


      #--------------------------------------------------------------------------------------------

      # CONDITIONAL SIMULATION AT GRIDDED LOCATIONS
      #--------------------------------------------

      else

        println(" ");
        println("Simulation at the grid nodes");
        println(" ");


        # How many rows can be simulated simultaneously?
        #-----------------------------------------------

        ntok = ceil(ntok/ng);
        m1 = max(1,min(ny,floor(ntok/nx)));
        rows = [λ for λ=0:m1:(ny-0.5)] ;
        push!(rows,ny)
        m2 = length(rows)-1;
        m3 = m2*nz;
        seed_nugget = abs.(Int.(rand(rng,Int64,m3)));


        # Loop over the grid rows
        #------------------------

        for n = 1:m3

          nnz = Int32(floor((n-1)/m2));
          nny = n-nnz*m2;
          index = [λ for λ=rows[nny]:(rows[nny+1]-1)];
          m4 = length(index); # number of rows to simulate
          m5 = m4*nx; # number of blocks to simulate
          m6 = m5*ng; # number of points discretizing the blocks


          # Coordinates of the points to simulate
          #--------------------------------------
          # println("ones(m4,1)*x0, y0 + index*dy ones(m4,1)*(z0 .+ nnnz*dz)=ones(",m4,",1),*",x0," ",y0," + ",index,"*",dy," ","ones(",m4,",1)*(",z0," .+ ",nnnz,"*",dz,")");
          coord0 = Int.(hcat(ones(m4,1)*x0, y0 + index*dy, ones(m4,1)*(z0 + nnz*dz)));
          #println("coord0=",coord0);
          # println("n=",n);
          # println("m4=",m4)
          blk_coord = kron(coord0,ones(nx,1)) + kron(ones(m4,1),Float64[λ for λ=0:(nx-1)]*[dx 0 0]); # coordinates of the block centres
          coord = kron(blk_coord,ones(ng,1))-kron(ones(m5,1),grid);
          # println("mean(coord)=",mean(coord))
          push!(COORDS,coord)

          # Non-conditional simulation
          #---------------------------

          simu = tbmain(coord,model,sigma,nlines,nrealiz,seed_line,all_lines,all_offset,
                        all_r,all_phi,all_theta,valid_lines);
          # println("G:mean(simu)=",mean(simu,1))

          # Add nugget effect
          #------------------

          if max_nugget > eps()
            # randn('state',seed_nugget[n]);
            rng = Random.MersenneTwister(seed_nugget[n]);
            simu = simu + (ones(m6,1)*total_nugget).*randn(rng,m6,nrealiz);
          end
          # println("H:mean(simu)=",mean(simu,1))
          # exit()
          # Conditioning
          #-------------

          if isinf(radius[1]) # dual kriging

            for i = 1:m6

              # Substitution of residuals
              k0 = setdual(model,coord[i,:],cc,b,datacoord,model_rotationmatrix,maxorder);
              simu[i,:] = simu[i,:] + k0*weights;

            end


          elseif radius[1] > eps() # kriging in a moving neighborhood

            for i = 1:m5

              # Coordinates of the centre of the ith block to simulate
              block_centre = [mean(coord[(i-1)*ng+1:i*ng,1]), mean(coord[(i-1)*ng+1:i*ng,2]), mean(coord[(i-1)*ng+1:i*ng,3])];
              # println(i);
              # println(mean(coord[(i-1)*ng+1:i*ng,1]), mean(coord[(i-1)*ng+1:i*ng,2]), mean(coord[(i-1)*ng+1:i*ng,3]))
              # # Search for neighboring data
              # println((#datacoord,residuals,
              # block_centre,search_rotationmatrix,octant,ndata,nxsup,
              #                  nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb));
              # println("mean(datacoord)=",mean(datacoord,1))
              # println("mean(residuals)=",mean(residuals,1))
              # println("mean(block_centre)=",mean(block_centre,1))
              # println("mean(search_rotationmatrix)=",mean(search_rotationmatrix,1))
              # println("octant=",octant)
              # println("ndata=",ndata)
              # println("nxsup=",nxsup)
              # println("nysup=",nysup)
              # println("nzsup=",nzsup)
              # println("xmnsup=",xmnsup)
              # println("ymnsup=",ymnsup)
              # println("zmnsup=",zmnsup)
              # println("xsizsup=",xsizsup)
              # println("ysizsup=",ysizsup)
              # println("zsizsup=",zsizsup)
              # println("mean(ixsbtosr)=",mean(ixsbtosr))
              # println("mean(iysbtosr)=",mean(iysbtosr))
              # println("mean(izsbtosr)=",mean(izsbtosr))
              # println("mean(nisb)=",mean(nisb))
              # println("size(nisb)=",size(nisb))
              # println("(nisb)=",(nisb))
              datacoord_i,residuals_i = Search(datacoord,residuals,block_centre,search_rotationmatrix,octant,ndata,nxsup,
                               nysup,nzsup,xmnsup,ymnsup,zmnsup,xsizsup,ysizsup,zsizsup,ixsbtosr,iysbtosr,izsbtosr,nisb);
              # Substitution of residuals
              if !isempty(datacoord_i) # samples are found in the kriging neighborhood
                weights = krige(datacoord_i,coord[(i-1)*ng+1:i*ng,:],model,cc,b,nugget+1e-7,model_rotationmatrix,maxorder);
                simu[(i-1)*ng+1:i*ng,:] = simu[(i-1)*ng+1:i*ng,:] + weights'*residuals_i;
                push!(RR,weights'*residuals_i);
              elseif (flag > 1) # intrinsic model, no conditioning data found
                println(" ");
                println("WARNING - No data found in moving neighborhood. Corresponding results may be senseless");
                println(" ");
                flag = flag - 1;
              end

            end

          end


          # Back transform to original scale
          #---------------------------------
          # println("tableZY=",tableZY)
          # println("I:mean(simu)=",mean(simu,1))
          # exit();
          simu = backtr(simu,tableZY,zmin,zmax,tail);
          # println("J:mean(simu)=",mean(simu,1))


          # Average the values over the blocks
          #-----------------------------------
          # println("J:size(simu)=",size(simu))

          if ng > 1
            # println(ng);
            simu = reshape(simu,(ng,m5*nrealiz));
            simu = mean(simu,1);
            simu = reshape(simu,(m5,Int64(nrealiz)));
          end
          # println("K:mean(simu)=",mean(simu,1))


          # Report on progress from time to time
          #-------------------------------------

          push!(progress,10*floor((10*n)/m3));
          if (progress[n+1] > progress[n])
            # disp(["  ",num2str(progress[n+1]),"# completed"]);
            sleep(1);
          end


          # Write in output file
          #---------------------

          # fprintf(fid,[outputformat,"\n"],simu');
          if n==1
              outputformat=simu
          else
              outputformat=vcat(outputformat,simu);
          end

        end

      end
      # println("F:size(simu)=",size(simu))
      return return_a_mapping_to_grades ? Float64.(outputformat) : α*log.(Float64.(outputformat))#,datacoord,extreme_coord,ydata,simucoord
      #--------------------------------------------------------------------------------------------

      # CLOSE THE OUTPUT FILE
      #----------------------

      # fclose(fid);
    end
  end

  export function tbsim_1(α=2.5,simucoord=Float64[])
    return tbsim(α,simucoord)
    # α=2.5,simucoord=Float64[],datacoord=[],ydata=[],nrealiz = 10,return_a_mapping_to_grades=true,nargin=2,nx=200,ny=300,nz=1,
    #     x0=1.0,y0=1.0,z0=100.0,dx=2.0,dy=2.0,dz=10.0,nd = Int16[1,1,1],limits=[-90,90],
    #     zmin=0.0,zmax = 10.0,tail=[1.0,5.0],nst=2,nugget=0.1,model=[],
    #     model1=[1 0.45 100 100 150 0 0 0 1 1000; 2 0.45 100 100 1000000000 0 0 0 1 1000],
    #     seed = 9784498, radius = [100 100 150] , angles =[0 0 0],
    #     octant = 1,ndata = 4,nbdecimal=3,ntok = 5000
    #
  end
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
