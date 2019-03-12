#----------------------------------------------------------------------------------
# Authors: Xavier Emery and Christian Lantuéjoul
# Paper: "TBSIM: A computer program for conditional simulation of three-dimensional
#         Gaussian random fields via the turning bands method"
#----------------------------------------------------------------------------------

function create_paramfile()
#--------------------------------------------------
# Create a default parameter file for program tbsim
#--------------------------------------------------
    open("tbsim.par", "w") do f
        write(f,"                  Parameters for TBSIM")
        write(f,"\n                   ********************")
        write(f,"\n  ")
        write(f,"\nSTART OF PARAMETERS:")
        write(f,"\n0                                % type of simulation: 0=gridded locations; 1=scattered locations")
        write(f,"\nlocations.prn                    % if =1: file with coordinates of locations to simulate")
        write(f,"\n1 2 3                            %        columns for location coordinates")
        write(f,"\n5.0 5.0 60.0                     % if =0: x0, y0, z0")
        write(f,"\n40 60 1                          %        nx, ny, nz")
        write(f,"\n10.0 10.0 12.0                   %        dx, dy, dz")
        write(f,"\n5 5 1                            % block discretization [1 1 1 for point-support simulation]")
        write(f,"\nnscore.out                       % file with conditioning data")
        write(f,"\n1 2 3                            %        columns for coordinates")
        write(f,"\n6                                %        column for Gaussian data")
        write(f,"\n-90  90                          %        trimming limits for Gaussian data")
        write(f,"\nnscore.trn                       % file with conversion table [raw-Gaussian]")
        write(f,"\n0.0 10.0                         % minimum and maximum values for untransformed variable")
        write(f,"\n1.0 5.0                          % parameters for lower()-tail and upper()-tail extrapolation")
        write(f,"\n2 0.1                            % number of nested structures, nugget effect")
        write(f,"\n2 0.45 100 100 150 0 0 0 1 1000  %        1st structure: it cc a1 a2 a3 ang1 ang2 ang3 b nlines")
        write(f,"\n2 0.45 100 100 1e9 0 0 0 1 1000  %        2nd structure: it cc a1 a2 a3 ang1 ang2 ang3 b nlines")
        write(f,"\n10                               % number of realizations")
        write(f,"\n9784498                          % seed for random number generation")
        write(f,"\n100 100 150                      % maximum search radii in the rotated system()")
        write(f,"\n0 0 0                            % angles for search ellipsoid()")
        write(f,"\n1                                % divide in octants? 1=yes, 0=no")
        write(f,"\n4                                % optimal number of data per octant [if octant=1] or in total [if 0]")
        write(f,"\ntbsim.out                        % name of output file")
        write(f,"\n3                                % number of decimals for values in the output file")
        write(f,"\n1                                % create a header in the output file? 1=yes, 0=no")
        write(f,"\n5000                             % maximum number of locations to simulate simultaneously")
        write(f,"\n ")
        write(f,"\nAvailable model types:")
        write(f,"\n        1: spherical")
        write(f,"\n        2: exponential")
        write(f,"\n        3: gamma (parameter b > 0)")
        write(f,"\n        4: stable [parameter b < 2]")
        write(f,"\n        5: cubic")
        write(f,"\n        6: Gaussian")
        write(f,"\n        7: cardinal sine")
        write(f,"\n        8: J-Bessel [parameter b > 0.5]")
        write(f,"\n        9: K-Bessel [parameter b > 0]")
        write(f,"\n       10: generalized Cauchy [parameter b > 0]")
        write(f,"\n       11: exponential sine")
        write(f,"\n       12: linear")
        write(f,"\n       13: power [exponent b > 0]")
        write(f,"\n       14: mixed power [exponent b <= 2]")
        write(f,"\n       15: spline (exponent b = even integer)")
    end
end
