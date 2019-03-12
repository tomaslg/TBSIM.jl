using PyPlot

#--------------------------------------------------------------------------------------------------------
# NON-CONDITIONAL SIMULATIONS
# TESTING THE REPRODUCTION OF THE VARIOGRAM
#------------------------------------------
#
# In this exercice; a set of 100 non-conditional realizations is constructed over a grid with size 50*50
# then their regional variograms are compared to the theoretical model.
#--------------------------------------------------------------------------------------------------------

clear()

# Grid parameters
simucoord = []
x0 = 0
y0 = 0
z0 = 0
nx = 150
ny = 150
nz = 1
dx = 1
dy = 1
dz = 1
nd = [1 1 1]; # point-support simulation
    
# We do not need data nor transformation table
datacoord = []
ydata = []
limits = [-90 90]
tableZY = []
zmin = 0
zmax = 1
tail = [1 1]

# Model parameters
it = 1;  # model type()
a1 = 10; # scale factor along y direction
a2 = 10; # scale factor along x direction
a3 = 10; # scale factor along z direction
ang1 = 0
ang2 = 0
ang3 = 0
model = [it a1 a2 a3 ang1 ang2 ang3]
cc = 0.9
b = .5
nugget = 0.1

# Simulation parameters
nlines = 1000
nrealiz = 100
seed = 9784498
ntok = 200

# We do not need neighborhood [non-conditional realizations]
radius = [0 0 0]
angles = [0 0 0]
octant = 0
ndata = 0

# Output file
name = "tbsim_nc.out"
nbdecimal = 3
header = 0

# Run the turning bands algorithm
tbsim(simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail, ...
      model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok)

# Check the simulated variograms along y direction
load() "tbsim_nc.out"; 
maxh = 20
clf()
for i = 1:nrealiz
  values = tbsim_nc[:,i]
  values = reshape(values,nx,ny)
  for h = 1:maxh
    D = 0.5*(values[:,1:ny-h]-values[:,1+h:ny]).^2
    variogram[i,h] = mean(D[:]); 
  end
  set(gcf(),"DefaultAxesFontName','Times','DefaultAxesFontSize",14)
  whitebg('w')
  hold on; 
  plot([1:maxh],variogram[i,:],"g--"); 
end
plot([1:maxh],mean(variogram),"r--")

# Superimpose the theoretical variogram model
C = cova[it,[1:maxh]/a1,b]
C0 = cova[it,0,b]
plot([1:maxh],nugget+cc*(C0-C),"b-")

pause()

#--------------------------------------------------------------------------------------------------------
# CONDITIONAL SIMULATIONS OF GRADES IN A BENCH OF A COPPER DEPOSIT
#-----------------------------------------------------------------
#
# This is an undocumented case study to test the post-processing of the realizations:
# conditioning to Gaussian data, back-transformation and change of support(). 
#
# The original data and normal scores transforms are available in the file nscore.out.
# The transformation table is contained in the file nscore.trn. Both files have been obtained thanks 
# to the GSLIB program nscore.exe. The normal scores variogram is a nested exponential + nugget model.
# The grid to simulate contains 40*60 blocks along the east and north directions, with a size of 10m*10m
# and a single block along the elevation with height 12m [length of the composite data].
#
#--------------------------------------------------------------------------------------------------------

clear()

# Grid parameters
simucoord = []
x0 = 5
y0 = 5
z0 = 60
nx = 40
ny = 60
nz = 1
dx = 10
dy = 10
dz = 12
nd = [5 5 1]

# Data and transformation table
load() "nscore.out"
datacoord = nscore[:,1:3]
ydata = nscore[:,6]
limits = [-90 90]
load() "nscore.trn"
tableZY = nscore
zmin = 0
zmax = 10
tail = [1 5]

# Model parameters
model = [2 100 100 150 0 0 0;2 100 100 1e5 0 0 0]
cc = [0.45;0.45]
b = [1;1]
nugget = 0.1

# Simulation parameters
nlines = 1000
nrealiz = 10
seed = 9784498
ntok = 500

# Neighborhood parameters
radius = [100 100 150]
angles = [0 0 0]
octant = 1
ndata = 4

# Output file
name = "tbsim_cu.out"
nbdecimal = 3
header = 0

# Run the turning bands algorithm
tbsim(simucoord,x0,y0,z0,nx,ny,nz,dx,dy,dz,nd,datacoord,ydata,limits,tableZY,zmin,zmax,tail, ...
      model,cc,b,nugget,nlines,nrealiz,seed,radius,angles,octant,ndata,name,nbdecimal,header,ntok)

# Display the realizations
load("tbsim_cu.out")
for i = 1:nrealiz
  values = tbsim_cu[:,i]
  values = reshape(values,nx,ny)
  figure(i)
  clf()
  set(gcf(),"DefaultAxesFontName','Times','DefaultAxesFontSize",14)
  whitebg('w')
  hold on; 
  pcolor([5:10:400],[5:10:600],min(3,values')); 
  colormap(jet)
  axis("image()")
  shading("flat")
  colorbar("vert")
end

#--------------------------------------------------------------------------------------------------------