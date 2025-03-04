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
fs = 22
font = {'family': 'Times New Roman',
                        'size'   : fs}

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

invRo = cdf_file.variables['R'][0]
invFr = np.sqrt(cdf_file.variables['B_therm'][0])
Fr = np.round(1/invFr,4)
Ro = np.round(1/invRo,4)

#these arrays are indexed by [t,z,y,x]
ux = np.array(cdf_file.variables['ux'][:])
uy = np.array(cdf_file.variables['uy'][:])
uz = np.array(cdf_file.variables['uz'][:])
temp = np.array(cdf_file.variables['Temp'][:])
wz = db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
uz_rms = np.sqrt(np.sum(uz**2,axis=1)/Nz)
wz_bar = np.sum(wz, axis=1)/Nz
tflux_bar = np.sum(uz*temp, axis=1)/Nz
wTmax = np.max(tflux_bar)
uz_rms_max = np.max(uz_rms)


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
alz = np.zeros([Nt, Nlz, Ny, Nx])

for i in range(Nlz):
    l = i
    lz[i] = dz*l
    alz[:,i,:,:] = ACFlz(uz, l)
    print("Completed Lengthscale: ", lz[i])

# perform FWHM
clz = np.zeros_like(alz[:,0,:,:])
for i in range(Nt):
    for j in range(Ny):
        for k in range(Nx):
            idx = np.where(alz[i,:,j,k] <= 0.5)
            if (len(idx[0]) == 0):
                clz[i,j,k] = lz[np.where(alz[i,:,j,k] == min(alz[i,:,j,k]))[0][0]]
            else:
                clz[i,j,k] = lz[idx[0][0]]

np.savez("vertavg", x = x, y = y, z = z, t = t, ux = ux,\
        uy = uy, wz = wz, wz_bar=wz_bar, tfluz_bar=tflux_bar, lz=lz, alz=alz,
        uz_rms=uz_rms)

fig, ax = plt.subplots(2, 2,figsize=(10,10))

# horizontal slider bar
#taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
#staxis = Slider(taxis, 'Height', 0, t[-1]-t[0], valinit=0)

# plot vertically averaged quantities
pc1 = ax[0,0].imshow(Fr*np.abs(wz_bar[0,:,:].T+invRo),
    norm=colors.Normalize(vmin=0,vmax=2), cmap='RdYlBu_r', origin='lower')
fig.colorbar(pc1, ax=ax[0,0])
ax[0,0].set_title(r"$Fr|\hat{\omega}_z+Ro^{-1}|$", **font)

pc2 = ax[0,1].imshow(tflux_bar[0,:,:].T,
    norm=colors.Normalize(vmin=-wTmax,vmax=wTmax), cmap='seismic', origin='lower')
fig.colorbar(pc2, ax=ax[0,1])
ax[0,1].set_title(r"$\hat{wT}$", **font)


pc3 = ax[1,0].imshow(clz[0,:,:].T,
    norm=colors.Normalize(vmin=dz,vmax=dz*Nlz), cmap='viridis', origin='lower')
fig.colorbar(pc3, ax=ax[1,0])
ax[1,0].set_title(r"$l_z$", **font)

pc4 = ax[1,1].imshow(clz[0,:,:].T,
    norm=colors.Normalize(vmin=0,vmax=uz_rms_max), cmap='plasma', origin='lower')
fig.colorbar(pc4, ax=ax[1,1])
ax[1,1].set_title(r"$\hat{u}_z$", **font)

for axisset in ax:
    for axis in axisset:
        axis.set_yticks([])
        axis.set_xticks([])

fig.suptitle("t = "+str(t[0])+", Fr = "+str(Fr)+", Ro = "+str(Ro), **font)
fig.tight_layout()

def update_frame(frame):
    pc1.set_array(Fr*np.abs(wz_bar[frame,:,:].T+invRo))
    pc2.set_array(tflux_bar[frame,:,:].T)
    pc3.set_array(clz[frame,:,:].T)
    pc4.set_array(uz_rms[frame,:,:].T)
    fig.suptitle("t = "+str(t[frame])+\
                ", Fr = "+str(Fr)+", Ro = "+str(Ro), **font)
    #staxis.set_val(t[frame]-t[0]) 
    print('Done with frame: ', frame)
    return (pc1, pc2, pc3, pc4)

ani = animation.FuncAnimation(fig=fig,
func=update_frame,frames=Nt,interval=1000,blit=True)
ani.save('RC_evolution.mp4')
ani.save('RC_evolution.gif')
