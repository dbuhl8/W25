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
#x = np.array(cdf_file.variables['x'])
#y = np.array(cdf_file.variables['y'])
#z = np.array(cdf_file.variables['z'])
#t = np.array(cdf_file.variables['t'][:])
#Nx = len(x)
#Ny = len(y)
#Nz = len(z)
#Nt = len(t)
#gx = cdf_file.variables['Gammax'][0]
#gy = cdf_file.variables['Gammay'][0]
#gz = cdf_file.variables['Gammaz'][0]
#dx = gx/Nx
#dy = gy/Ny
#dz = gz/Nz

invRo = cdf_file.variables['R']
invFr = np.sqrt(cdf_file.variables['B_therm'])

#these arrays are indexed by [t,z,y,x]
#ux = np.array(cdf_file.variables['ux'][:])
#uy = np.array(cdf_file.variables['uy'][:])
#uz = np.array(cdf_file.variables['uz'][:])
#temp = np.array(cdf_file.variables['Temp'][:])
#wz =  db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
#tdisp = np.sqrt(FD6X(temp,Nx,dx)**2 +\
        #FD6Y(temp,Ny,dy)**2 + \
        #FD4Z(temp,Nz,dz)**2)
#uxmax = np.max(np.abs(ux),axis=0)
#uymax = np.max(np.abs(uy),axis=0)
#uzmax = np.max(np.abs(uz),axis=0)
#wzmax = np.max(np.abs(wz),axis=0)/4

#tempmax = np.max(np.abs(temp),axis=0)
#tdispmax = np.max(np.abs(tdisp),axis=0)

# plotting with matplotlib
#XX, YY = np.meshgrid(x,y) # if spatial discretization is needed
#fig, ax = plt.subplots(2, 2)

# compute vertical averages
#ux_bar = np.sum(ux,axis=1)/Nz
#uy_bar = np.sum(uy,axis=1)/Nz
#uz_bar = np.sum(uz,axis=1)/Nz   
#temp_bar = np.sum(temp,axis=1)/Nz 
#wz_bar = np.sum(wz,axis=1)/Nz   

#np.savez("Om8B100_vertavg", x = x, y = y, z = z, t = t, ux = ux,\
        #uy = uy, wz = wz, wz_bar=wz_bar)
npz_file = np.load("vertavg.npz")
t = npz_file['t']
wz_bar = npz_file['wz_bar']
tflux_bar= npz_file['tflux_bar']
ux = npz_file['ux']
uy = npz_file['uy']
wTmax = np.max(tflux_bar)
Nt = len(t)
lz = npz_file['lz']
Nlz = len(lz)
alz = npz_file['alz']

# compute FWHM

fig, ax = plt.subplots(1, 2,figsize=(16,9))

# horizontal slider bar
taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
staxis = Slider(taxis, 'Height', 0, t[-1]-t[0], valinit=0)

# plot vertically averaged quantities
pc1 = ax[0].imshow(invFr/np.maximum(wz_bar[0,:,:].T+invRo,1e-5),
    norm=colors.Normalize(vmin=0,vmax=1), cmap='plasma_r', origin='lower')
fig.colorbar(pc1, ax=ax[0])
ax[0].set_title(r"$\frac{Fr^{-1}}{\hat{\omega}_z+Ro^{-1}}$")

pc2 = ax[1].imshow(tflux_bar[0,:,:].T,
    norm=colors.Normalize(vmin=-wTmax,vmax=wTmax), cmap='RdYlBu_r', origin='lower')
fig.colorbar(pc2, ax=ax[1])
ax[1].set_title(r"$\hat{\omega}_z+Ro^{-1}$")


def update_frame(frame):
    pc1.set_array(invFr/np.maximum(wz_bar[frame,:,:].T+invRo, 1e-5))
    pc2.set_array(tflux_bar[frame,:,:].T)
    staxis.set_val(t[frame]-t[0]) 
    print('Done with frame: ', frame)
    return (pc1, pc2)

ani = animation.FuncAnimation(fig=fig,
func=update_frame,frames=Nt,interval=1000,blit=True)
ani.save('RC_evolution.gif')

#plt.show()

