import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import xarray as xr

# this file is meant to work for a simdat.cdf file. 
# different configurations exsit for slice files which are 2D in nature. 
def FDY(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,:,0] = (field[:,:,1]-field[:,:,0])/dy
    dy_field[:,:,ny-1] = (field[:,:,ny-1]-field[:,:,ny-2])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,:,i+1] = (0.5/dy)*(field[:,:,i+2]-field[:,:,i])
    return dy_field

def FDZ(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dz_field[:,0,:] = (field[:,1,:]-field[:,0,:])/dz
    dz_field[:,nz-1,:] = (field[:,nz-1,:]-field[:,nz-2,:])/dz
    # centered finite difference formula (inside)
    for i in range(nz-2):
        dz_field[:,i+1,:] = (0.5/dz)*(field[:,i+2,:]-field[:,i,:])
    return dz_field

fn = 'YZSLICE1.cdf'
cdf_file = xr.open_dataset(fn)
dtype = np.float32

#obtaining discretization data
y = cdf_file.y.to_numpy().astype(dtype)
z = cdf_file.z.to_numpy().astype(dtype)
t = cdf_file.t.to_numpy().astype(dtype)
Ny = len(y)
Nz = len(z)
Nt = len(t)
gy = 4*np.pi
gz = np.pi
dy = gy/Ny
dz = gz/Nz
dt = cdf_file.dt.to_numpy().astype(dtype)

#these arrays are indexed by [t,z,y,x]
ux = cdf_file.ux.to_numpy().astype(dtype)
uy = cdf_file.uy.to_numpy().astype(dtype)
uz = cdf_file.uz.to_numpy().astype(dtype)
temp = cdf_file.Temp.to_numpy().astype(dtype)
wz = FDY(uy, Ny, dy) - FDZ(ux, Nz, dz)
tdisp = np.sqrt(FDY(temp, Ny, dy)**2 \
                + FDZ(temp, Nz, dz)**2)
wzmax = 50
uxmax = np.max(np.abs(ux[:,:,:]))
uymax = np.max(np.abs(uy[:,:,:]))
uzmax = np.max(np.abs(uz[:,:,:]))
tempmax = np.max(np.abs(temp[:,:,:]))

# plotting with matplotlib
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

# plot using imshow (much faster than pcolor for evenly spaced discretization
# axis style is worse however
pc1 = ax[0,0].imshow(ux[0,:,:],
    norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc1, ax=ax[0,0], orientation='horizontal')
ax[0,0].set_title(r'$u_x$')

pc2 = ax[0,1].imshow(uy[0,:,:],
    norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic',
    origin='lower')
fig.colorbar(pc2, ax=ax[0,1], orientation='horizontal')
ax[0,1].set_title(r'$u_y$')

pc3 = ax[1,0].imshow(temp[0,:,:],
    norm=colors.Normalize(vmin=-tempmax,vmax=tempmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc3, ax=ax[1,0], orientation='horizontal')
ax[1,0].set_title(r"$T'$")

pc4 = ax[1,1].imshow(uz[0,:,:],
    norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic',
    origin='lower')
fig.colorbar(pc4, ax=ax[1,1], orientation='horizontal')
ax[1,1].set_title(r'$u_z$')

# makes room for the slider
fig.subplots_adjust(0.05, .1, .95, .9,
.2,.2)

# horizaontally oriented slider
taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
staxis = Slider(taxis, 'Time', t[0], t[-1],
valinit=t[0], valstep=dt)

# vertically orientated slider
#zaxis = plt.axes([0.94,0.15,0.03,0.7],facecolor='blue',orientation='vertical')
#szaxis = Slider(zaxis, 'Height', 0, np.pi, valinit=0, valstep=dz)

for axis_set in ax:
    for axis in axis_set:
        axis.set_xticks([])
        axis.set_yticks([])

#fig.tight_layout()

# movie down vertical extent of domain
def update_frame(frame):
    pc1.set_array(ux[frame,:,:])
    pc2.set_array(uy[frame,:,:])
    pc3.set_array(temp[frame,:,:])
    pc4.set_array(uz[frame,:,:])
    staxis.set_val(t[frame]) 
    print('Done with frame: ', frame)
    return (pc1,pc2,pc3,pc4)

def slide_update(val):
    return update_frame(np.argmax(t == val))

ani = animation.FuncAnimation(fig=fig,
func=update_frame,frames=Nt,interval=100,blit=True)
ani.save('YZ_evolution.gif')
#plt.show()

