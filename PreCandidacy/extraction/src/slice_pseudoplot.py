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
    # returns dataset (numpy) of same shape but computing the 1st derivative # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,0,:] = (field[:,1,:]-field[:,0,:])/dy
    dy_field[:,ny-1,:] = (field[:,ny-1,:]-field[:,ny-2,:])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,i+1,:] = (0.5/dy)*(field[:,i+2,:]-field[:,i,:])
    return dy_field

def make_norm(maxnum):
    return colors.Normalize(vmin=-maxnum,vmax=maxnum)

def sclmap(cmap, maxnum):
    return mpl.cm.ScalarMappable(cmap=cmap,norm=make_norm(maxnum))

# opening slice data files
fnxy = 'XYSLICE1.cdf'
fnxz = 'XZSLICE1.cdf'
fnyz = 'YZSLICE1.cdf'
cdf_filexy = xr.open_dataset(fnxy)
cdf_filexz = xr.open_dataset(fnxz)
cdf_fileyz = xr.open_dataset(fnyz)
dtype = np.float32
cmap = 'RdYlBu_r'
num_contours = 100 # enforces smoothness

# obtaining discretization data
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

uy_xy = cdf_filexy.uy.to_numpy().astype(dtype) # [t, y, x]
uy_xz = cdf_filexz.uy.to_numpy().astype(dtype) # [t, z, x]
uy_yz = cdf_fileyz.uy.to_numpy().astype(dtype) # [t, z, y]

temp_xy = cdf_filexy.Temp.to_numpy().astype(dtype) # [t, y, x]
temp_xz = cdf_filexz.Temp.to_numpy().astype(dtype) # [t, z, x]
temp_yz = cdf_fileyz.Temp.to_numpy().astype(dtype) # [t, z, y]

uz_xy = cdf_filexy.uz.to_numpy().astype(dtype) # [t, y, x]
uz_xz = cdf_filexz.uz.to_numpy().astype(dtype) # [t, z, x]
uz_yz = cdf_fileyz.uz.to_numpy().astype(dtype) # [t, z, y]

# Finding max values for the colorbar
uxmax = np.max([ux_xz.max(), ux_yz.max(), ux_xy.max()])
uymax = np.max([uy_xz.max(), uy_yz.max(), ux_xy.max()])
uzmax = np.max([uz_xz.max(), uz_yz.max(), uz_xy.max()])
tempmax = np.max([temp_xz.max(), temp_yz.max(), temp_xy.max()])

#print('Maximum; ux: ', uxmax, ', uy: ',\
    #uymax, ', uz: ', uzmax, ', temp: ',tempmax)

kwux = {
    'vmin': -uxmax,
    'vmax': uxmax,
    #'norm': make_norm(uxmax),
    'cmap': cmap
}

kwuy = {
    'vmin': -uymax,
    'vmax': uymax,
    #'norm': make_norm(uymax),
    'cmap': cmap
}

kwuz = {
    'vmin': -uzmax,
    'vmax': uzmax,
    #'norm': make_norm(tempmax),
    'cmap': cmap
}

kwtemp = {
    'vmin': -tempmax,
    'vmax': tempmax,
    #'norm': make_norm(uzmax),
    'cmap': cmap
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
    ------------------------------------------------------------------------------  
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
    -----------------------------------------------------------------------------
    need to implement this for each subplot
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
fig, ax = plt.subplots(2, 2,subplot_kw=dict(projection='3d'))

# ux box plot
ux_top = ax[0,0].contourf(XX[:,:,0],YY[:,:,0],ux_xy[0,:,:], levels=num_contours,
    zdir='z', offset=0, **kwux)
ux_xzside = ax[0,0].contourf(XX[0,:,:],ux_xz[0,:,:].T,
    ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwux)
ux_yzside = ax[0,0].contourf(ux_yz[0,:,:].T,YY[:,-1,:],
    ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwux)
ax[0,0].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
ax[0,0].set_title(r"$u_x$")
ax[0,0].view_init(40, 240, 0)
ax[0,0].set_box_aspect((4, 4, 1), zoom=0.9)
fig.colorbar(sclmap(cmap,uxmax), ax=ax[0,0])

# uy box plot
uy_top = ax[0,1].contourf(XX[:,:,0],YY[:,:,0],uy_xy[0,:,:], levels=num_contours,
    zdir='z', offset=0, **kwuy)
uy_xzside = ax[0,1].contourf(XX[0,:,:],uy_xz[0,:,:].T,
    ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuy)
uy_yzside = ax[0,1].contourf(uy_yz[0,:,:].T,YY[:,-1,:],
    ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuy)
ax[0,1].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
ax[0,1].set_title(r"$u_y$")
ax[0,1].view_init(40, 240, 0)
ax[0,1].set_box_aspect((4, 4, 1), zoom=0.9)
#fig.colorbar(uy_top, ax=ax[0,1], norm=make_norm(uymax), fraction=0.02, pad=0.1)
fig.colorbar(sclmap(cmap,uymax), ax=ax[0,1])


# temp box plot
temp_top = ax[1,0].contourf(XX[:,:,0],YY[:,:,0],temp_xy[0,:,:],
    levels=num_contours, zdir='z', offset=0, **kwtemp)
temp_xzside = ax[1,0].contourf(XX[0,:,:],temp_xz[0,:,:].T,
    ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwtemp)
temp_yzside = ax[1,0].contourf(temp_yz[0,:,:].T,YY[:,-1,:],
    ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwtemp)
ax[1,0].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
ax[1,0].set_title(r"$T'$")
ax[1,0].view_init(40, 240, 0)
ax[1,0].set_box_aspect((4, 4, 1), zoom=0.9)
#fig.colorbar(temp_top, ax=ax[1,0], norm=make_norm(tempmax), fraction=0.02, pad=0.1)
fig.colorbar(sclmap(cmap,tempmax), ax=ax[1,0])


# uz box plot
uz_top = ax[1,1].contourf(XX[:,:,0],YY[:,:,0],uz_xy[0,:,:], levels=num_contours,
    zdir='z', offset=0, **kwuz)
uz_xzside = ax[1,1].contourf(XX[0,:,:],uz_xz[0,:,:].T,
    ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuz)
uz_yzside = ax[1,1].contourf(uz_yz[0,:,:].T,YY[:,-1,:],
    ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuz)
ax[1,1].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
ax[1,1].set_title(r"$u_z$")
ax[1,1].view_init(40, 240, 0)
ax[1,1].set_box_aspect((4, 4, 1), zoom=0.9)
#fig.colorbar(uz_top, ax=ax[1,1], norm=make_norm(uzmax), fraction=0.02, pad=0.1)
fig.colorbar(sclmap(cmap,uzmax), ax=ax[1,1])


# Additional Plot Formating
# makes room for the slider
fig.subplots_adjust(0, .05, .90, .90,
.05,.05)
for axis_set in ax:
    for axis in axis_set:
        axis.set_xticks([])
        axis.set_yticks([])
        axis.set_zticks([])
# horizaontally oriented slider
taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
staxis = Slider(taxis, 'Time', t[0], t[-1], valinit=t[0], valstep=dt)

# movie down vertical extent of domain
def update_frame(frame):
    ux_top = ax[0,0].contourf(XX[:,:,0],YY[:,:,0],
        ux_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwux)
    ux_xzside = ax[0,0].contourf(XX[0,:,:],ux_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwux)
    ux_yzside = ax[0,0].contourf(ux_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwux)
    uy_top = ax[0,1].contourf(XX[:,:,0],YY[:,:,0],
        uy_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwuy)
    uy_xzside = ax[0,1].contourf(XX[0,:,:],uy_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuy)
    uy_yzside = ax[0,1].contourf(uy_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuy)
    temp_top = ax[1,0].contourf(XX[:,:,0],YY[:,:,0],
        temp_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwtemp)
    temp_xzside = ax[1,0].contourf(XX[0,:,:],temp_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwtemp)
    temp_yzside = ax[1,0].contourf(temp_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwtemp)
    uz_top = ax[1,1].contourf(XX[:,:,0],YY[:,:,0],
        uz_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwuz)
    uz_xzside = ax[1,1].contourf(XX[0,:,:],uz_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuz)
    uz_yzside = ax[1,1].contourf(uz_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuz)
    staxis.set_val(t[frame]) 
    print('Done with frame: ', frame)
    return (ux_top, ux_xzside, ux_yzside, uy_top, uy_xzside, uy_yzside,
        temp_top, temp_xzside, temp_yzside, uz_top, uz_xzside, uz_yzside)

ani = animation.FuncAnimation(fig=fig,
    func=update_frame,frames=Nt,interval=100,blit=True)
ani.save('psuedoplot.gif')
#plt.show()

