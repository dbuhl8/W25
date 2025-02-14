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

fn = 'simdat1.cdf'
cdf_file = xr.open_dataset(fn)
dtype = np.float32

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
#dX = np.zeros(3)
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz
#print('Gridstep: ',dx, dy, dz)

#these arrays are indexed by [t,z,y,x]
ux = cdf_file.ux.to_numpy().astype(dtype)
uy = cdf_file.uy.to_numpy().astype(dtype)
uz = cdf_file.uz.to_numpy().astype(dtype)
wz = FDX(uy, Nx, dx) - FDY(ux, Ny, dy)
wzmax = 50 #np.max(np.abs(wz[0,0,:,:]))
uxmax = np.max(np.abs(ux[-1,:,:,:]))
uymax = np.max(np.abs(uy[-1,:,:,:]))
uzmax = np.max(np.abs(uz[-1,:,:,:]))

# plotting with matplotlib
XX, YY = np.meshgrid(x,y)
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
pc1 = ax[0,0].imshow(ux[-1,0,:,:].T,
    norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc1, ax=ax[0,0])
ax[0,0].set_title(r'$u_x$')

pc2 = ax[0,1].imshow(uy[-1,0,:,:].T,
    norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic',
    origin='lower')
fig.colorbar(pc2, ax=ax[0,1])
ax[0,1].set_title(r'$u_y$')

pc3 = ax[1,0].imshow(wz[-1,0,:,:].T,
    norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc3, ax=ax[1,0])
ax[1,0].set_title(r'$\omega_z$')

pc4 = ax[1,1].imshow(uz[-1,0,:,:].T,
    norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc4, ax=ax[1,1])
ax[1,1].set_title(r'$u_z$')

fig.subplots_adjust(0.1, .1, .9, .9,
.1,.1)

zaxis = plt.axes([0.1, 0.02, 0.9, 0.03], facecolor='blue')
szaxis = Slider(zaxis, 'Height', 0, np.pi,
valinit=0, valstep=dz)

#fig.tight_layout()

# movie down vertical extent of domain
def update_frame(frame):
    #uxmax = np.max(np.abs(ux[-1,frame,:,:]))
    #uymax = np.max(np.abs(uy[-1,frame,:,:]))
    #uzmax = np.max(np.abs(uz[-1,frame,:,:]))
    pc1.set_array(ux[-1,frame,:,:].T)
    pc2.set_array(uy[-1,frame,:,:].T)
    pc3.set_array(wz[-1,frame,:,:].T)
    pc4.set_array(uz[-1,frame,:,:].T)
    szaxis.set_val(frame*dz) 
    #pc1.set_norm(colors.Normalize(vmin=-uxmax,vmax=uxmax))
    #pc2.set_norm(colors.Normalize(vmin=-uymax,vmax=uymax))
    #pc4.set_norm(colors.Normalize(vmin=-uzmax,vmax=uzmax))
    print('Done with frame: ', frame)
    return (pc1,pc2,pc3,pc4)

#def update_frame(frame):
    #uxmax = np.max(np.abs(ux[-1,frame,:,:]))
    #uymax = np.max(np.abs(uy[-1,frame,:,:]))
    #uzmax = np.max(np.abs(uz[-1,frame,:,:]))
    #pc1.set_array(ux[-1,frame,:,:].T)
    #pc2.set_array(uy[-1,frame,:,:].T)
    #pc3.set_array(wz[-1,frame,:,:].T)
    #pc4.set_array(uz[-1,frame,:,:].T)
    #pc1.set_norm(colors.Normalize(vmin=-uxmax,vmax=uxmax))
    #pc2.set_norm(colors.Normalize(vmin=-uymax,vmax=uymax))
    #pc4.set_norm(colors.Normalize(vmin=-uzmax,vmax=uzmax))
    #return [[pc1,pc2],[pc3,pc4]]

ani = animation.FuncAnimation(fig=fig,
func=update_frame,frames=Nz,interval=50,blit=True)
ani.save('imshow_z_extent.gif')

