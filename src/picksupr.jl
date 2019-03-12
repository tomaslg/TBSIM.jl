#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function picksupr(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,search_rotationmatrix)

#-------------------------------------------------------------------------------------
# Build template of super-blocks centered at the block containing the node to simulate
#-------------------------------------------------------------------------------------

# Main Loop over all possible super blocks

      nsbtosr = 0;
      ixsbtosr=Int32[];
      iysbtosr=Int32[];
      izsbtosr=Int32[];
      for i=-nxsup+1:nxsup-1

        for j=-nysup+1:nysup-1

          for k = -nzsup+1:nzsup-1

            xo = i*xsizsup;
            yo = j*ysizsup;
            zo = k*zsizsup;
            coord = abs.(Float64[xo yo zo]);

            # Calculate the closest distance between the corners of the super block and the block at the origin

            distance = [coord[1]-xsizsup, coord[2]-ysizsup, coord[3]-zsizsup]';
            distance = distance*search_rotationmatrix;
            distance = sum(distance'.^2)';

            # Keep this super block if it is close enough:

            if (distance < 1)
              nsbtosr = nsbtosr + 1;
              push!(ixsbtosr,i); # ixsbtosr[nsbtosr] = i;
              push!(iysbtosr,j); # iysbtosr[nsbtosr] = j;
              push!(izsbtosr,k); # izsbtosr[nsbtosr] = k;
            end

          end

        end

      end
return ixsbtosr,iysbtosr,izsbtosr
end
