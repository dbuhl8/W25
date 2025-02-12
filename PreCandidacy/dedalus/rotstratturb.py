
import numpy as np
import dedalus.public as d3
import logging

logger = logging.getLogger(__name__)


# Parameters 
Lx, Ly, Lz = 4*np.pi, 4*np.pi, np.pi
Nx, Ny, Nz = 256, 256, 64
Re, Pe, Ro, Fr = 600, 60, 1./3, 0.1
dealias = 3/2
dtype = np.float64
timestepper = d3.timesteppers.RK443
max_timestep = 1e-2

# Bases (what is a distributor?)
coords = d3.CartesianCoordinates('x', 'y', 'z')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias = dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias = dealias)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias = dealias)

# Fields 
# (looks like the distributor defines variables over the coordinate bases)
p = dist.Field(name='p', bases=(xbasis, ybasis, zbasis))
t = dist.Field(name='t', bases=(xbasis, ybasis, zbasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis,zbasis))
f = dist.VectorField(coords, name='f', bases=(xbasis,ybasis,zbasis))
# need to transcribe gaussian regression code to here

# Substitutions
nu = 1/Re
kappa = 1/Pe
x, y, z = dist.local_grids(xbasis,ybasis,zbasis)
ex, ey, ez = coords.unit_vector_fields(dist)

# Problem
#problem = d3.IVP([u, t, p, f,

# need to finish writing the equations (need to ensure periodic BC)
# need to complete forcing mechanism
# need to do some other stuff as well
