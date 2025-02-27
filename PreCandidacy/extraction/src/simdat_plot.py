import numpy as np
import dbuhlMod as db
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from netCDF4 import MFDataset

fn = 'simdat*.cdf'
cdf_file = MFDataset(fn)
dtype = np.float32
#ptstp = -1  # this choses which timestep to plot (-1 is the last one)

#obtaining discretization data
x = np.array(cdf_file.variables['x'])
y = np.array(cdf_file.variables['y'])
z = np.array(cdf_file.variables['z'])
t = np.array(cdf_file.variables['t'][:])
Nx = len(x)
Ny = len(y)
Nz = len(z)
Nt = len(t)
gx = cdf_file.variables['Gammax'][0]
gy = cdf_file.variables['Gammay'][0]
gz = cdf_file.variables['Gammaz'][0]
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz

invRo = cdf_file.variables['R']
invFr = np.sqrt(cdf_file.variables['B_therm'])

#these arrays are indexed by [t,z,y,x]
ux = np.array(cdf_file.variables['ux'][:])
uy = np.array(cdf_file.variables['uy'][:])
uz = np.array(cdf_file.variables['uz'][:])
wz = db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
wz_bar = np.sum(wz, axis=1)/Nz


def ACFlz(field, l):
    # performs integral of f(z)f(z+l) over vertical domain
    # and normalizes it by f^2(z)
    num_z = len(field[0,:,0,0])
    norm = np.sum(field**2, axis=1)
    val = np.zeros_like(norm)
    for fj in range(num_z):
        val += field[:,(fj+l)%num_z,:,:]*field[:,fj,:,:]
    norm = np.divide(val,np.maximum(norm, 1e-5))
    return norm


Nlz = int(Nz/2)
lz = np.zeros(Nlz)
idx = np.where(invFr/(wz_bar+invRo) < 1)
alz = np.zeros([Nt, Nlz, Ny, Nz])

for i in range(Nlz):
    l = Nlz - i
    lz[i] = dz*l
    alz[:,i,:,:] = ACFlz(uz, l)
    print("Completed Lengthscale: ", lz[i])

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
"""

np.savez("vertavg", x = x, y = y, z = z, t = t, ux = ux,\
        uy = uy, wz = wz, wz_bar=wz_bar, lz=lz, alz=alz)

