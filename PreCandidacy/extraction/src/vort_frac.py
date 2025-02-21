import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from netCDF4 import Dataset
#import xarray as xr

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


def FD4X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,1] = (field[:,:,:,2]-field[:,:,:,0])/(2*dx)
    dx_field[:,:,:,nx-2] = (field[:,:,:,nx-1]-field[:,:,:,nx-3])/(2*dx)
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # 4th order centered finite difference formula (inside)
    for i in range(nx-4):
        dx_field[:,:,:,i+2] = (1/(12*dx))*(-field[:,:,:,i+4]+8*field[:,:,:,i+3]\
                                           -8*field[:,:,:,i+1]+field[:,:,:,i])
    return dx_field

def FD4Y(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # 4th order centered finite difference formula (inside)
    for i in range(ny-4):
        dy_field[:,:,i+2,:] = (1/(12*dy))*(-field[:,:,i+4,:]+8*field[:,:,i+3,:]\
                                           -8*field[:,:,i+1,:]+field[:,:,i,:])
    return dy_field

def FD4Z(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    dz_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dz
    dz_field[:,1,:,:] = (field[:,2,:,:]-field[:,0,:,:])/(2*dz)
    dz_field[:,nz-2,:,:] = (field[:,nz-1,:,:]-field[:,nz-3,:,:])/(2*dz)
    dz_field[:,nz-1,:,:] = (field[:,nz-1,:,:]-field[:,nz-2,:,:])/dz
    # 4th order centered finite difference formula (inside)
    for i in range(nz-4):
        dz_field[:,i+2,:,:] = (1/(12*dz))*(-field[:,i+4,:,:]+8*field[:,i+3,:,:]\
                                           -8*field[:,i+1,:,:]+field[:,i,:,:])
    return dz_field

def FD6X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,1] = (field[:,:,:,2]-field[:,:,:,0])/(2*dx)
    dx_field[:,:,:,2] = (-field[:,:,:,4]+8*field[:,:,:,3] \
                        -8*field[:,:,:,1]+field[:,:,:,0])/(12*dx)
    dx_field[:,:,:,nx-3] = (-field[:,:,:,nx-1]+8*field[:,:,:,nx-2]\
                            -8*field[:,:,:,nx-5]+field[:,:,:,nx-6])/(12*dx)
    dx_field[:,:,:,nx-2] = (field[:,:,:,nx-1]-field[:,:,:,nx-3])/(2*dx)
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # 4th order centered finite difference formula (inside)
    for i in range(nx-6):
        dx_field[:,:,:,i+3] = (1/(60*dx))*(-field[:,:,:,i+6]+9*field[:,:,:,i+5]\
                                           -45*field[:,:,:,i+4]\
                                           +45*field[:,:,:,i+2]-9*field[:,:,:,i+1]\
                                           +field[:,:,:,i])
    return -dx_field

def FD6Y(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    dy_field[:,:,2,:] = (-field[:,:,4,:]+8*field[:,:,3,:] \
                        -8*field[:,:,1,:]+field[:,:,0,:])/(12*dy)
    dy_field[:,:,ny-3,:] = (-field[:,:,ny-1,:]+8*field[:,:,ny-2,:]\
                            -8*field[:,:,ny-5,:]+field[:,:,ny-6,:])/(12*dy)
    dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # 4th order centered finite difference formula (inside)
    for i in range(ny-6):
        dy_field[:,:,i+3,:] = (1/(60*dy))*(-field[:,:,i+6,:]+9*field[:,:,i+5,:]\
                                           -45*field[:,:,i+4,:]\
                                           +45*field[:,:,i+2,:]-9*field[:,:,i+1,:]\
                                           +field[:,:,i,:])
    return -dy_field

# commented out because this breaks expanse plotting
#plt.rcParams['text.usetex'] = True

dtype = np.float32
ptstp = -1  # this choses which timestep to plot (-1 is the last one)
fs = 18
dfs = 9
font = {'family': 'Times New Roman',
                        'size'   : fs}

files_30 = ['Om0.5B30Re600Pe60/OUT.dat',\
    'Om1B30Re600Pe60/OUT.dat','Om2B30Re600Pe60/OUT.dat', \
    'Om3B30Re600Pe60/OUT.dat','Om4B30Re600Pe60/OUT.dat',\
    'Om5B30Re600Pe60/OUT.dat','Om8B30Re600Pe60/OUT.dat','Om10B30Re600Pe60/OUT.dat']

files_100 = ['Om0.5B100Re600Pe60/OUT.dat',\
    'Om1B100Re600Pe60/OUT.dat','Om2B100Re600Pe60/OUT.dat', \
    'Om3B100Re600Pe60/OUT.dat','Om4B100Re600Pe60/OUT.dat',\
    'Om8B100Re600Pe60/OUT.dat']
cor_idx = [0, 1, 2, 3, 4, 6]

simdat_files_30 = ['Om0.5B30Re600Pe60/simdat4.cdf',\
    'Om1B30Re600Pe60/simdat3.cdf','Om2B30Re600Pe60/simdat4.cdf',\
    'Om3B30Re600Pe60/simdat7.cdf','Om4B30Re600Pe60/simdat4.cdf',\
    'Om5B30Re600Pe60/simdat2.cdf', 'Om8B30Re600Pe60/simdat1.cdf',\
    'Om10B30Re600Pe60/simdat4.cdf']
sim_30_idx = [0, 1, 2, 3, 4, 5, 6, 7]

simdat_files_100 = ['Om0.5B100Re600Pe60/simdat3.cdf',\
    'Om1B100Re600Pe60/simdat2.cdf','Om2B100Re600Pe60/simdat4.cdf', \
    'Om3B100Re600Pe60/simdat4.cdf','Om4B100Re600Pe60/simdat2.cdf',\
    'Om8B100Re600Pe60/simdat5.cdf']
sim_100_idx = [0,1,2,3,4,5]

invRo_30 = [0.5, 1, 2, 3, 4, 5, 8, 10]
invRo_100 = [0.5, 1, 2, 3, 4, 8]
invFr = [np.sqrt(30), 10]

labels_100 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=8$']
labels_30 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=5$',r'$\Omega=8$',r'$\Omega=10$']
cor_label_30 = [1,2,3,4,5,7]


num_files = [len(files_30), len(files_100)]
sim_num_files = [len(simdat_files_30), len(simdat_files_100)]
num_samples = 5 # num of samples from each timeseries
cvar = 0.2

vol_frac_30 = np.zeros(sim_num_files[0])
avg_uh_30 = np.zeros_like(vol_frac_30)
vol_frac_100 = np.zeros(sim_num_files[1])
avg_uh_100 = np.zeros_like(vol_frac_100)

for i, fn in enumerate(simdat_files_30):
    #cdf_file = xr.open_dataset(fn)
    cdf_file = Dataset(fn)

    #obtaining discretization data
    x = np.array(cdf_file.variables["x"])
    y = np.array(cdf_file.variables["y"])
    z = np.array(cdf_file.variables["z"])
    Nx = len(x)
    Ny = len(y)
    Nz = len(z)
    gx = cdf_file.variables["Gammax"][0]
    gy = cdf_file.variables["Gammay"][0]
    gz = cdf_file.variables["Gammaz"][0]
    dx = gx/Nx
    dy = gy/Ny
    dz = gz/Nz

    #these arrays are indexed by [t,z,y,x]
    ux = np.array(cdf_file.variables["ux"])
    uy = np.array(cdf_file.variables["uy"])
    avg_uh_30[i] = np.sqrt(np.sum(ux[ptstp,:,:,:]**2 + uy[ptstp,:,:,:]**2))
    wz =  FD6X(uy, Nx, dx) - FD6Y(ux, Ny, dy)
    wzmax = np.max(np.abs(wz[ptstp,:,:,:]))

    #compute volume fraction where wz >= .1 wzmax
    vol_frac_30[i] = len(np.where(np.abs(invFr[0]/(wz[ptstp,:,:,:]+invRo_30[i]))\
                        < 1)[0])/len(wz[ptstp,:,:,:])
    print("Finished with file: ", fn, ".")
    cdf_file.close()

for i, fn in enumerate(simdat_files_100):
    cdf_file = Dataset(fn)

    x = np.array(cdf_file.variables["x"])
    y = np.array(cdf_file.variables["y"])
    z = np.array(cdf_file.variables["z"])
    Nx = len(x)
    Ny = len(y)
    Nz = len(z)
    gx = cdf_file.variables["Gammax"][0]
    gy = cdf_file.variables["Gammay"][0]
    gz = cdf_file.variables["Gammaz"][0]
    dx = gx/Nx
    dy = gy/Ny
    dz = gz/Nz

    #these arrays are indexed by [t,z,y,x]
    ux = np.array(cdf_file.variables["ux"])
    uy = np.array(cdf_file.variables["uy"])


    avg_uh_100[i] = np.sqrt(np.sum(ux[ptstp,:,:,:]**2 + uy[ptstp,:,:,:]**2))
    wz =  FD6X(uy, Nx, dx) - FD6Y(ux, Ny, dy)
    wzmax = np.max(np.abs(wz[ptstp,:,:,:]))/4

    #compute volume fraction where wz >= .1 wzmax
    vol_frac_100[i] = len(np.where(\
                np.abs(invFr[1]/(wz[ptstp,:,:,:]+invRo_100[i]))\
                        < 1)[0])/len(wz[ptstp,:,:,:])
    print("Finished with file: ", fn, ".")
    cdf_file.close()

eInvRo_30 = np.zeros([num_files[0], num_samples])
eInvRo_err_30 = np.zeros_like(eInvRo_30)
avg_wT_30 = np.zeros_like(eInvRo_30)
avg_wT_err_30 = np.zeros_like(eInvRo_30)

eInvRo_100 = np.zeros([num_files[1], num_samples])
eInvRo_err_100 = np.zeros_like(eInvRo_100)
avg_wT_100 = np.zeros_like(eInvRo_100)
avg_wT_err_100 = np.zeros_like(eInvRo_100)

i = 0
for fn in files_30:
    file = open(fn, 'r')
    urms = np.array([])
    t = np.array([])
    wT = np.array([])
    # load data from text file to np arrays
    for line in file:
        if not line.startswith('#') and line.strip():  # skips headers
            # save data from this line (can choose which data to save)
            line_elements = [float(x) for x in line.split()]
            t = np.append(t,line_elements[1])
            urms = np.append(urms,line_elements[3])
            wT = np.append(wT,line_elements[8])

    # perform averaging over eddie turnover times (roughly)
    dt = 2*(np.pi**2)/urms[0]
    idx = np.nonzero(t == t[0])[0]
    for j in range(num_samples):
        cent = t[idx[-1]] + dt/2
        idx = np.nonzero(np.abs(t-cent) <= dt/2)[0]
        eInvRo_30[i,j] = np.sum(1/urms[idx])*invRo_30[i]/len(idx)
        eInvRo_err_30[i,j] = np.std(urms[idx])
        avg_wT_30[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err_30[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
    i += 1
    print("Finished with file: ", fn, ".")
    file.close()

i = 0
for fn in files_100:
    file = open(fn, 'r')
    urms = np.array([])
    t = np.array([])
    wT = np.array([])
    # load data from text file to np arrays
    for line in file:
        if not line.startswith('#') and line.strip():  # skips headers
            # save data from this line (can choose which data to save)
            line_elements = [float(x) for x in line.split()]
            t = np.append(t,line_elements[1])
            urms = np.append(urms,line_elements[3])
            wT = np.append(wT,line_elements[8])

    # perform averaging over eddie turnover times (roughly)
    dt = 2*(np.pi**2)/urms[0]
    idx = np.nonzero(t == t[0])[0]
    for j in range(num_samples):
        cent = t[idx[-1]] + dt/2
        idx = np.nonzero(np.abs(t-cent) <= dt/2)[0]
        eInvRo_100[i,j] = np.sum(1/urms[idx])*invRo_100[i]/len(idx)
        eInvRo_err_100[i,j] = np.std(urms[idx])
        avg_wT_100[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err_100[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
    i += 1
    print("Finished with file: ", fn, ".")
    file.close()

fig, ax = plt.subplots(2, 2,figsize=(16,9))
# rearrange ax arary to be indexed 1:4
#ax = [axitem for axitem in [axset for axset in ax]] 
dc = cvar/num_files[0]
col = np.linspace(0+dc,1-dc,num_files[0], True)
colors = [np.linspace(c-dc,c+dc,num_samples, True) for c in col]

artists = []

for i in range(num_files[0]):
    for j in range(num_samples):
        item = ax[0,0].errorbar(eInvRo_30[i,j], avg_wT_30[i,j],
        avg_wT_err_30[i,j], eInvRo_err_30[i,j],\
        'none', cm.nipy_spectral_r(colors[i][j]))
        if j==2:
            artists.append(item)

ax[0,0].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[0,0].set_ylabel(r'$wT$',rotation=0, **font)
ax[0,0].tick_params(axis='both', labelsize=fs-dfs)
ax[0,0].set_yscale('log')
ax[0,0].set_xscale('log')
ax[0,0].set_title(r'$Fr = 0.18$', **font)
ax[0,0].legend(artists,labels_30,loc='lower left',fontsize=fs-dfs)

#dc = cvar/num_files[1]
#col = np.linspace(0+dc,1-dc,num_files[1], True)
#colors = [np.linspace(c-dc,c+dc,num_samples, True) for c in col]
artists = []

for i in range(num_files[1]):
    for j in range(num_samples):
        item = ax[0,1].errorbar(eInvRo_100[i,j], avg_wT_100[i,j],
        avg_wT_err_100[i,j], eInvRo_err_100[i,j],\
        'none', cm.nipy_spectral_r(colors[cor_idx[i]][j]))
        if j==2:
            artists.append(item)

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
ax[0,1].set_xlabel(r'$Ro_e^{-1}$', **font)
#ax[1].set_ylabel(r'$wT$',rotation=0,labelpad=30, **font)
#ax[1].set_ylabel(r'$wT$',rotation=0, **font)
ax[0,1].tick_params(axis='both', labelsize=fs-dfs)
ax[0,1].sharey(ax[0,0])
ax[0,1].sharex(ax[0,0])
ax[0,1].set_yscale('log')
ax[0,1].set_xscale('log')
ax[0,1].set_title(r'$Fr = 0.1$', **font)
ax[0,1].legend(artists,labels_100,loc='lower left',fontsize=fs-dfs)

for i in range(sim_num_files[0]):
    artists.append(ax[1,0].scatter(invRo_30[i]/avg_uh_30[i], vol_frac_30[i], 
    c=cm.nipy_spectral_r(colors[sim_30_idx[i]][2])))

ax[1,0].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[1,0].set_ylabel(r'$\frac{1}{|V|}\int_{Fr/Ro_e < 1}dV$',rotation=0, **font)
ax[1,0].tick_params(axis='both', labelsize=fs-dfs)
ax[1,0].set_yscale('log')
ax[1,0].set_xscale('log')
ax[1,0].set_title(r'$Fr = 0.18$', **font)
ax[1,0].legend(artists,labels_30,loc='lower left',fontsize=fs-dfs)

artists = []

for i in range(sim_num_files[1]):
    artists.append(ax[1,1].scatter(invRo_100[i]/avg_uh_100[i],
    vol_frac_100[i], 
    c=cm.nipy_spectral_r(colors[cor_idx[i]][2])))

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
ax[1,1].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[1,1].tick_params(axis='both', labelsize=fs-dfs)
ax[1,1].sharey(ax[1,0])
ax[1,1].sharex(ax[1,0])
ax[1,1].set_yscale('log')
ax[1,1].set_xscale('log')
ax[1,1].set_title(r'$Fr = 0.1$', **font)
ax[1,1].legend(artists,labels_100,loc='lower left',fontsize=fs-dfs)

#fig.tight_layout()
#plt.show()
plt.savefig('flux_vol_frac_plot.pdf')

# plotting with matplotlib
#XX, YY = np.meshgrid(x,y) # if spatial discretization is needed
#fig, ax = plt.subplots(2, 2)

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

""" z extent gif plot 
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

