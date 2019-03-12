#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function tbmain(coord,model,sigma,nlines,nrealiz,seed,all_lines,all_offset,all_r,all_phi,all_theta,valid_lines)

#------------------------------------------------------------
# Non conditional simulation by the turning bands method
# Main routine for program tbsim (simulation along the lines)
#------------------------------------------------------------

# Define parameters
#------------------

  nst = size(model,1);
  m = size(coord,1);

# Initialize output
#------------------

  simu = zeros(m,nrealiz);

# Loop over each nested structure
#--------------------------------

  for i = 1:nst

    if (model[i,1] < 4.5) | (model[i,1] == 9)
        # Spherical, exponential, gamma, stable (b<1) and K-Bessel (b<0.5) models

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];
        valid = findall( xy-> xy > 0, valid_lines[index,i]);
        nbvalid = size(valid,1);

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        x = coord*lines';
        offset = ones(m,1)*all_offset[index,i]';
        interval_number = floor.(x-offset) .+ 1;
        position = x-offset - interval_number .+ 0.5;
        state = seed[(k-1)*nst+i,1:nlines[i]];
        # println("size(position)=",size(position))
        # println("size(valid)=",size(valid))
        # println("size(interval_number)=",size(interval_number))
        # Loop over the lines
        #--------------------

        for j = 1:nbvalid

          # Simulate the values by a partition method
          #------------------------------------------
          rng = Random.MersenneTwister(state[valid[j]]);
          maxnum = Int64(maximum(interval_number[:,valid[j]]));
          slope = sigma[i]*(2*(randn(rng,maxnum).>0) .- 1);
          simu[:,k] = simu[:,k] + slope[Int64.(interval_number[:,valid[j]])].*position[:,valid[j]];
          # println("size(valid[j])=",size(valid[j]))
          # println("size(position[:,valid[j]])=",size(position[:,valid[j]]))
          # println("size(Int64.(interval_number[:,valid[j]]))=",size(Int64.(interval_number[:,valid[j]])))
          # println("size(slope[Int64.(interval_number[:,valid[j]])].*position[:,valid[j]])=",size(slope[Int64.(interval_number[:,valid[j]])].*position[:,valid[j]]))
        end
        # exit()
      end

    elseif (model[i,1] == 5) # Cubic model

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];
        valid = Int.(findall( xy-> xy > 0, valid_lines[index,i]));
        nbvalid = size(valid,1);

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        x = coord*lines';
        offset = ones(m,1)*all_offset[index,i]';
        interval_number = Int.(floor.(x-offset) .+ 1);
        position = x-offset - interval_number .+ 0.5;
        state = seed[(k-1)*nst+i,1:nlines[i]];


        # Loop over the lines
        #--------------------

        for j = 1:nbvalid

          # Simulate the values by a partition method
          #------------------------------------------
          rng = Random.MersenneTwister(state[valid[j]]);
          maxnum = Int(maximum(interval_number[:,valid[j]]));
          slope = sigma[i]*(2*(randn(rng,maxnum).>0) .- 1);
          simu[:,k] = simu[:,k] + slope[interval_number[:,valid[j]]].*(0.25*position[:,valid[j]] - position[:,valid[j]].^3);

        end

      end

    elseif (model[i,1] < 12) # Stable (b>1), Gaussian, cardinal sine, Bessel J, Bessel K (b>0.5),
                             # generalized Cauchy and exponential sine models

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        lines = all_lines[((k-1)*nlines[i]+1):(k*nlines[i]),:,i];
        x = coord*lines';

        # Simulate the values by the continuous spectral method
        #------------------------------------------------------

        r = all_r[:,k,i] * ones(1,m);
        phi = all_phi[:,k,i] * ones(1,m);
        simu[:,k] = simu[:,k] + sigma[i]*sum(cos.(r.*x' .+ phi),dims=1)';

      end


    elseif (model[i,1] == 12) # Linear model

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        x = coord*lines';
        offset = ones(m,1)*all_offset[index,i]';
        interval_number = Int.(floor.(x-offset) .+ 1);
        state = seed[(k-1)*nst+i,1:nlines[i]];

        # Loop over the lines
        #--------------------

        for j = 1:nlines[i]

          # Simulate the values by a partition method
          #------------------------------------------

          rng = Random.MersenneTwister(state[j]);
          maxnum = maximum(interval_number[:,j]);
          stairs = sigma[i]*cumsum(2*(randn(rng,maxnum).>0) .- 1);
          simu[:,k] = simu[:,k] + stairs[interval_number[:,j]];

        end

      end


    elseif (model[i,1] == 13) # Power model (b<=1)

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        x = coord*lines';
        offset = ones(m,1)*all_offset[index,i]';
        interval_number = Int.(floor.(x-offset) .+ 1);
        state = seed[(k-1)*nst+i,1:nlines[i]];

        # Loop over the lines
        #--------------------

        for j = 1:nlines[i]

          # Simulate the values by a partition method
          #------------------------------------------
          rng = Random.MersenneTwister(state[j]);
          maxnum = length(interval_number[:,j]);#maximum(interval_number[:,j]);
          value = sigma[i]*randn(rng,maxnum);
          simu[:,k] = simu[:,k] + value;#[interval_number[:,j]];

        end

      end


    elseif (model[i,1] == 13.1) # Pure drift model

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];
        x = coord*lines';

        # Simulate the values by a monomial
        #----------------------------------

        x = x.^all_phi[1,1,i];
        r = all_r[:,k,i] * ones(1,m);
        simu[:,k] = simu[:,k] + sigma[i]*sum(r.*x',1)';

      end


    elseif (model[i,1] >= 13.2) # Power (b>1), mixed power and spline models

      # Loop over the realizations
      #---------------------------

      for k = 1:nrealiz

        # Project the points to simulate over the lines of the i-th nested structure
        #---------------------------------------------------------------------------

        index = ((k-1)*nlines[i]+1):(k*nlines[i]);
        lines = all_lines[index,:,i];
        x = coord*lines';

        # Simulate the values by a continuous spectral method
        #----------------------------------------------------

        r = all_r[:,k,i] * ones(1,m);
        phi = all_phi[:,k,i] * ones(1,m);
        theta = all_theta[:,k,i] * ones(1,m);
        simu[:,k] = simu[:,k] + sigma[i]* sum( theta.*cos.(r.*x' .+ phi), dims=1  )';

      end

    end

  end
  return simu
end
