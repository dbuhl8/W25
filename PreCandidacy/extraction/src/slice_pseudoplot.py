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
    dx_field[:,:,0] = (field[:,:,1]-field[:,:,0])/dx
    dx_field[:,:,nx-1] = (field[:,:,nx-1]-field[:,:,nx-2])/dx
    # centered finite difference formula (inside)
    for i in range(nx-2):
        dx_field[:,:,i+1] = (0.5/dx)*(field[:,:,i+2]-field[:,:,i])
    return dx_field

def FDY(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,0,:] = (field[:,1,:]-field[:,0,:])/dy
    dy_field[:,ny-1,:] = (field[:,ny-1,:]-field[:,ny-2,:])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,i+1,:] = (0.5/dy)*(field[:,i+2,:]-field[:,i,:])
    return dy_field

fnxy = 'XYSLICE1.cdf'
fnxz = 'XZSLICE1.cdf'
fnyz = 'YZSLICE1.cdf'
cdf_filexy = xr.open_dataset(fnxy)
cdf_filexz = xr.open_dataset(fnxz)
cdf_fileyz = xr.open_dataset(fnyz)
dtype = np.float32

#obtaining discretization data
x = cdf_filexy.x.to_numpy().astype(dtype)
y = cdf_filexy.y.to_numpy().astype(dtype)
z = cdf_filexz.z.to_numpy().astype(dtype)
t = cdf_filexy.t.to_numpy().astype(dtype)
Nx = len(x)
Ny = len(y)
Nz = len(z)
Nt = len(t)
gx = 4*np.pi
gy = 4*np.pi
gz = np.pi
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz
dt = cdf_filexy.dt.to_numpy().astype(dtype)

#these arrays are indexed by [t,z,y,x]
ux_xy = cdf_filexy.ux.to_numpy().astype(dtype) # [t, y, x]
ux_xz = cdf_filexz.ux.to_numpy().astype(dtype) # [t, z, x]
ux_yz = cdf_fileyz.ux.to_numpy().astype(dtype) # [t, z, y]
uy_xy = cdf_filexy.uy.to_numpy().astype(dtype)
uz_xy = cdf_filexy.uz.to_numpy().astype(dtype)
temp = cdf_filexy.Temp.to_numpy().astype(dtype)
#wz = FDX(uy, Nx, dx) - FDY(ux, Ny, dy)
#tdisp = np.sqrt(FDX(temp, Nx, dx)**2 \
                #+ FDY(temp, Ny, dy)**2)
#wzmax = 50
#uxmax = np.max(np.abs(ux[:,:,:]))
#uymax = np.max(np.abs(uy[:,:,:]))
#uzmax = np.max(np.abs(uz[:,:,:]))

uxmax = np.max([ux_xz.max(), ux_yz.max(), ux_xy.max()])

kw = {
    'vmin': -uxmax,
    'vmax': uxmax,
    'cmap': 'RdYlBu_r'
}

# plotting with matplotlib
XX, YY, ZZ = np.meshgrid(x,y, -z)

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
    # plot using imshow (much faster than pcolor for evenly spaced discretization)
    # axis style is worse however
    pc1 = ax[0,0].imshow(ux[0,:,:].T,
        norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc1, ax=ax[0,0])
    ax[0,0].set_title(r'$u_x$')

    pc2 = ax[0,1].imshow(uy[0,:,:].T,
        norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc2, ax=ax[0,1])
    ax[0,1].set_title(r'$u_y$')

    pc3 = ax[1,0].imshow(wz[0,:,:].T,
        norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc3, ax=ax[1,0])
    ax[1,0].set_title(r'$\omega_z$')

    pc4 = ax[1,1].imshow(uz[0,:,:].T,
        norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc4, ax=ax[1,1])
    ax[1,1].set_title(r'$u_z$')
"""

""" need to implement this for each subplot
    import matplotlib.pyplot as plt
    import numpy as np

    # Define dimensions
    Nx, Ny, Nz = 100, 300, 500
    X, Y, Z = np.meshgrid(np.arange(Nx), np.arange(Ny), -np.arange(Nz))

    # Create fake data
    data = (((X+100)**2 + (Y-20)**2 + 2*Z)/1000+1)

    kw = {
        'vmin': data.min(),
        'vmax': data.max(),
        'levels': np.linspace(data.min(), data.max(), 10),
    }

    # Create a figure with 3D ax
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111, projection='3d')

    # Plot contour surfaces
    _ = ax.contourf(
        X[:, :, 0], Y[:, :, 0], data[:, :, 0],
        zdir='z', offset=0, **kw
    )
    _ = ax.contourf(
        X[0, :, :], data[0, :, :], Z[0, :, :],
        zdir='y', offset=0, **kw
    )
    C = ax.contourf(
        data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
        zdir='x', offset=X.max(), **kw
    )
    # --


    # Set limits of the plot from coord limits
    xmin, xmax = X.min(), X.max()
    ymin, ymax = Y.min(), Y.max()
    zmin, zmax = Z.min(), Z.max()
    ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

    # Plot edges
    edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
    ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

    # Set labels and zticks
    ax.set(
        xlabel='X [km]',
        ylabel='Y [km]',
        zlabel='Z [m]',
        zticks=[0, -150, -300, -450],
    )

    # Set zoom and angle view
    ax.view_init(40, -30, 0)
    ax.set_box_aspect(None, zoom=0.9)

    # Colorbar
    fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Name [units]')

    # Show Figure
    plt.show()
"""

# makes a box plot for each subplot
#fig, ax = plt.subplots(2, 2, projection='3d')
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111, projection='3d')

# ux plot
# top of box (XY SLICE)
C = ax.contourf(XX[:,:,0],YY[:,:,0],
ux_xy[0,:,:],zdir='z', offset=0, **kw)
# x side of box (XZ SLICE)
ax.contourf(XX[0,:,:],ux_xz[0,:,:].T,
ZZ[0,:,:], zdir='y', offset=0, **kw)
# y-side of box (YZ SLICE)
ax.contourf(ux_yz[0,:,:].T,YY[:,-1,:],
ZZ[:,-1,:],zdir='x', offset=0, **kw)

ax.set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

#ax.view_init(40, -30, 0)
#ax.set_box_aspect(None, zoom=0.9)

fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1)
plt.show()

# makes room for the slider
#fig.subplots_adjust(0.05, .1, .95, .9,
#.2,.2)

# horizaontally oriented slider
#taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
#staxis = Slider(taxis, 'Time', t[0], t[-1],
#valinit=t[0], valstep=dt)

# vertically orientated slider
#zaxis = plt.axes([0.94,0.15,0.03,0.7],facecolor='blue',orientation='vertical')
#szaxis = Slider(zaxis, 'Height', 0, np.pi, valinit=0, valstep=dz)

#for axis_set in ax:
    #for axis in axis_set:
        #axis.set_xticks([])
        #axis.set_yticks([])

# movie down vertical extent of domain
def update_frame(frame):
    pc1.set_array(ux[frame,:,:].T)
    pc2.set_array(uy[frame,:,:].T)
    pc3.set_array(wz[frame,:,:].T)
    pc4.set_array(uz[frame,:,:].T)
    staxis.set_val(t[frame]) 
    print('Done with frame: ', frame)
    return (pc1,pc2,pc3,pc4)

#ani = animation.FuncAnimation(fig=fig,
#func=update_frame,frames=Nt,interval=100,blit=True)
#ani.save('XY_evolution.gif')

