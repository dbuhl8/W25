import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import xarray as xr

# this file is meant to work for a simdat.cdf file. 
# different configurations exsit for slice files which are 2D in nature. 
def FDX(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # centered finite difference formula (inside)
    for i in range(nx-2):
        dx_field[:,:,:,i+1] = (0.5/dx)*(field[:,:,:,i+2]-field[:,:,:,i])
    return dx_field

def FDY(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,:,i+1,:] = (0.5/dy)*(field[:,:,i+2,:]-field[:,:,i,:])
    return dy_field

def FDZ(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dz_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dz
    dz_field[:,nz-1,:,:] = (field[:,nz-1,:,:]-field[:,nz-2,:,:])/dz
    # centered finite difference formula (inside)
    for i in range(nz-2):
        dz_field[:,i+1,:,:] = (0.5/dz)*(field[:,i+2,:,:]-field[:,i,:,:])
    return dz_field

fn = 'simdat1.cdf'
cdf_file = xr.open_dataset(fn)
dtype = np.float32
ptstp = -1 # this choses which timestep to plot (-1 is the last one)

#obtaining discretization data
x = cdf_file.x.to_numpy()
y = cdf_file.y.to_numpy()
z = cdf_file.z.to_numpy()
Nx = len(x)
Ny = len(y)
Nz = len(z)
#print('Number of gridpoints: ',Nx,Ny,Nz)
gx = cdf_file.Gammax.to_numpy().astype(dtype)
gy = cdf_file.Gammay.to_numpy().astype(dtype)
gz = cdf_file.Gammaz.to_numpy().astype(dtype)
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz
#print('Gridstep: ',dx, dy, dz)

#these arrays are indexed by [t,z,y,x]
ux = cdf_file.ux.to_numpy().astype(dtype)
uy = cdf_file.uy.to_numpy().astype(dtype)
uz = cdf_file.uz.to_numpy().astype(dtype)
temp = cdf_file.Temp.to_numpy().astype(dtype)
wz = FDX(uy, Nx, dx) - FDY(ux, Ny, dy)
tdisp = np.sqrt(FDX(temp,Nx,dx)**2 + FDY(temp,Ny,dy)**2 + FDZ(temp,Nz,dz)**2)
uxmax = np.max(np.abs(ux[ptstp,:,:,:]))
uymax = np.max(np.abs(uy[ptstp,:,:,:]))
uzmax = np.max(np.abs(uz[ptstp,:,:,:]))
wzmax = 2*max([uxmax, uymax, uzmax]) #np.max(np.abs(wz[0,0,:,:]))
tempmax = np.max(np.abs(temp[ptstp,:,:,:]))
tdispmax = np.max(np.abs(tdisp[ptstp,:,:,:]))

# plotting with matplotlib
#XX, YY = np.meshgrid(x,y) # if spatial discretization is needed
fig, ax = plt.subplots(2, 2)

"""
    # can loop over this to create a movie through the vertical domain
    pc1 = ax[0,0].pcolor(XX[:,:],
        YY[:,:], ux[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic')
    fig.colorbar(pc1, ax=ax[0,0])
    ax[0,0].set_title(r'$u_x$')

    pc2 = ax[0,1].pcolor(XX[:,:],
        YY[:,:], uy[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic')
    fig.colorbar(pc2, ax=ax[0,1])
    ax[0,1].set_title(r'$u_y$')

    pc3 = ax[1,0].pcolor(XX[:,:],
        YY[:,:], wz[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic')
    fig.colorbar(pc3, ax=ax[1,0])
    ax[1,0].set_title(r'$\omega_z$')

    pc4 = ax[1,1].pcolor(XX[:,:],
        YY[:,:], uz[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic')
    fig.colorbar(pc4, ax=ax[1,1])
    ax[1,1].set_title(r'$u_z$')
"""
# can loop over this to create a movie through the vertical domain
pc1 = ax[0,0].imshow(temp[ptstp,0,:,:].T,
    norm=colors.Normalize(vmin=-tempmax,vmax=tempmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc1, ax=ax[0,0])
ax[0,0].set_title(r"$T'$")

pc2 = ax[0,1].imshow(tdisp[ptstp,0,:,:].T,
    norm=colors.Normalize(vmin=0,vmax=tdispmax), cmap='viridis',
    origin='lower')
fig.colorbar(pc2, ax=ax[0,1])
ax[0,1].set_title(r'$|\nabla T|^2$')

pc3 = ax[1,0].imshow(wz[ptstp,0,:,:].T,
    norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc3, ax=ax[1,0])
ax[1,0].set_title(r'$\omega_z$')

pc4 = ax[1,1].imshow(uz[ptstp,0,:,:].T,
    norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc4, ax=ax[1,1])
ax[1,1].set_title(r'$u_z$')

# adjusts spacing so that slider fits
fig.subplots_adjust(0.01, .1, .91, .9,
.2,.15)

# horizontal slider bar
#zaxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
#szaxis = Slider(zaxis, 'Height', 0, np.pi,
#valinit=0, valstep=dz)

# vertical slider bar
zaxis = plt.axes([0.93, 0.15, 0.03, 0.7], facecolor='blue')
szaxis = Slider(zaxis, 'Height', 0, np.pi,
valinit=0, valstep=dz, orientation='vertical')

for axis_set in ax:
    for axis in axis_set:
        axis.set_xticks([])
        axis.set_yticks([])


# movie down vertical extent of domain
def update_frame(frame):
    pc1.set_array(temp[ptstp,frame,:,:].T)
    pc2.set_array(tdisp[ptstp,frame,:,:].T)
    pc3.set_array(wz[ptstp,frame,:,:].T)
    pc4.set_array(uz[ptstp,frame,:,:].T)
    szaxis.set_val(frame*dz) 
    print('Done with frame: ', frame)
    return (pc1,pc2,pc3,pc4)

ani = animation.FuncAnimation(fig=fig,
func=update_frame,frames=Nz,interval=50,blit=True)
ani.save('imshow_z_extent.gif')

