from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from netCDF4 import Dataset

fn = "XZSLICE1.cdf"

# open cdf file as a dataset
cdf_file = Dataset(fn)

# indexed by [t,z,x]
ux = np.array(cdf_file.variables["ux"])
uy = np.array(cdf_file.variables["uy"])
uz = np.array(cdf_file.variables["uz"])
temp = np.array(cdf_file.variables["Temp"])

x = np.array(cdf_file.variables["x"])
z = np.array(cdf_file.variables["z"])
t = np.array(cdf_file.variables["t"])

Nx = len(x)
Nz = len(z)
Nt = len(t)
dx = 4*np.pi/Nx
dz = np.pi/Nz
dt = np.array(cdf_file.variables["timestep"])

bux = np.sum(ux, axis=1)/Nz
buy = np.sum(uy, axis=1)/Nz

def ACFlz(field, l):
    # performs integral of f(z)f(z+l) over vertical domain
    # and normalizes it by f^2(z)
    num_t = len(field[:,0,0])
    num_z = len(field[0,:,0])
    num_x = len(field[0 ,0,:])
    norm = np.zeros([num_t,num_x])
    norm = np.sum(field**2, axis=1)
    val = np.zeros_like(norm)
    for fj in range(num_z):
        val += field[:,(fj+l)%num_z,:]*field[:,fj,:]
    norm = np.divide(val,np.maximum(norm, 1e-5))
    return norm

Nlz = int(Nz/4)

lz= np.zeros(Nlz)
alz= np.zeros([4, Nt, Nlz, Nx])


for i in range(len(z)):
    ux[:,i,:] -= bux
    uy[:,i,:] -= buy

for i in range(Nlz):
    l = Nlz - i
    lz[i] = dz*l
    alz[0,:,i,:] = ACFlz(ux,l)
    alz[1,:,i,:] = ACFlz(uy,l)
    alz[2,:,i,:] = ACFlz(uz,l)
    alz[3,:,i,:] = ACFlz(temp,l)
    print("Completed Lengthscale: ", lz[i])

np.savez("lz_data", lz = lz, alz = alz, ux = ux, uy = uy, x = x)

ptstp = 27
invRo = 2
Fr = np.round(1/np.sqrt(30),4)
Ro = np.round(1/invRo,4)


#npzfile = np.load("lz_data.npz")
#lz = npzfile["lz"]
#alz = npzfile["alz"]
#ux = npzfile["ux"]
#uy = npzfile["uy"]
#x = npzfile["x"]



uxmax = np.max(np.abs(ux[ptstp,:,:]))
uymax = np.max(np.abs(uy[ptstp,:,:]))
umax = max([uxmax, uymax])

fig, ax = plt.subplots(3, 2)

im1 = ax[0,0].imshow(alz[0,ptstp,:,:],vmin=0,vmax=1,cmap="plasma",aspect='auto')
im2 = ax[0,1].imshow(alz[1,ptstp,:,:],vmin=0,vmax=1,cmap="plasma",aspect='auto')
im3 = ax[1,0].imshow(alz[2,ptstp,:,:],vmin=0,vmax=1,cmap="plasma",aspect='auto')
im4 = ax[1,1].imshow(alz[3,ptstp,:,:],vmin=0,vmax=1,cmap="plasma",aspect='auto')
im5 = ax[2,0].imshow(ux[ptstp,:,:],vmin=-umax,vmax=umax,cmap="RdYlBu_r")
im6 = ax[2,1].imshow(uy[ptstp,:,:],vmin=-umax,vmax=umax,cmap="RdYlBu_r")

ax[0,0].set_title(r"$ACF(u_x')$")
ax[0,1].set_title(r"$ACF(u_y')$")
ax[1,0].set_title(r'$ACF(u_z)$')
ax[1,1].set_title(r'$ACF(temp)$')
ax[2,0].set_title(r"$u_x'$")
ax[2,1].set_title(r"$u_y'$")
ax[2,0].set_ylabel(r'$Z$')
ax[2,1].set_ylabel(r'$Z$')

ntic = 8
f1 = int(Nlz/ntic)

ax[0,0].set_yticks(np.arange(f1, Nlz, f1))
ax[0,1].set_yticks(np.arange(f1, Nlz, f1))
ax[1,0].set_yticks(np.arange(f1, Nlz, f1))
ax[1,1].set_yticks(np.arange(f1, Nlz, f1))
ax[0,0].set_yticklabels(np.round(np.flip(np.arange(f1, Nlz, f1))*dz,2))
ax[0,1].set_yticklabels(np.round(np.flip(np.arange(f1, Nlz, f1))*dz,2))
ax[1,0].set_yticklabels(np.round(np.flip(np.arange(f1, Nlz, f1))*dz,2))
ax[1,1].set_yticklabels(np.round(np.flip(np.arange(f1, Nlz, f1))*dz,2))

ax[2,0].set_yticks([])
ax[2,1].set_yticks([])


fig.subplots_adjust(right=0.8,hspace=.4,wspace=.3)
cbar_ax = fig.add_axes([0.85, 0.40, 0.05, 0.45])
cbar_ax2 = fig.add_axes([0.85, 0.10, 0.05, 0.25])
fig.colorbar(im4, cax=cbar_ax)
fig.colorbar(im5, cax=cbar_ax2)
fig.suptitle(r't_i = '+str(t[ptstp])+', i = '+str(ptstp)+', Fr = '+str(Fr)+', Ro = '+str(Ro))

for axis_set in ax:
    for axis in axis_set:
        axis.set_xticks([])
        #axis.set_xlabel("X")

def update_frame(frame):
    im1.set_array(alz[0,frame,:,:])
    im2.set_array(alz[1,frame,:,:])
    im3.set_array(alz[2,frame,:,:])
    im4.set_array(alz[3,frame,:,:])
    im5.set_array(ux[frame,:,:])
    im6.set_array(uy[frame,:,:])
    fig.suptitle(r't_i = '+str(t[frame])+', i = '+str(frame)+', Fr = '+str(Fr))
    print('Done with frame: ', frame)
    return (im1,im2,im3,im4,im5,im6)

ani = animation.FuncAnimation(fig=fig,
    func=update_frame,frames=Nt,interval=500,blit=True)
ani.save('acflz_evolution.mp4')
#plt.show()
