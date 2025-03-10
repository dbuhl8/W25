import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.legend_handler import HandlerTuple
from matplotlib.widgets import Slider
from netCDF4 import Dataset

# commented out because this breaks expanse plotting
plt.rcParams['text.usetex'] = True

dtype = np.float32
ptstp = -1  # this choses which timestep to plot (-1 is the last one)
fs = 18
dfs = 6
ms = 20
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

simdat_files_30 = ['Om0.5B30Re600Pe60/simdat*.cdf',\
    'Om1B30Re600Pe60/simdat*.cdf','Om2B30Re600Pe60/simdat*.cdf',\
    'Om3B30Re600Pe60/simdat*.cdf','Om4B30Re600Pe60/simdat*.cdf',\
    'Om5B30Re600Pe60/simdat*.cdf', 'Om8B30Re600Pe60/simdat*.cdf',\
    'Om10B30Re600Pe60/simdat*.cdf']
sim_30_idx = [0, 1, 2, 3, 4, 5, 6, 7]

simdat_files_100 = ['Om0.5B100Re600Pe60/simdat*.cdf',\
    'Om1B100Re600Pe60/simdat*.cdf','Om2B100Re600Pe60/simdat*.cdf', \
    'Om3B100Re600Pe60/simdat*.cdf','Om4B100Re600Pe60/simdat*.cdf',\
    'Om8B100Re600Pe60/simdat*.cdf']
sim_100_idx = [0,1,2,3,4,5]

invRo_30 = [0.5, 1, 2, 3, 4, 5, 8, 10]
invRo_100 = [0.5, 1, 2, 3, 4, 8]
invFr = [np.sqrt(30), 10]

labels_100 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=8$']
labels_30 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=5$',r'$\Omega=8$',r'$\Omega=10$']
cor_label_30 = [1,2,3,4,5,7]


num_files = [len(files_30), len(files_100)]
sim_num_files = [len(simdat_files_30), len(simdat_files_100)]
cvar = 0.3

# there is a section missing from vort_frac.py where data extraction is
# performed 

# this is the command used to save the data
#np.savez("vort_frac", eInvRo_30=eInvRo_30, avg_wT_30=avg_wT_30,\
    #eInvRo_err_30=eInvRo_err_30, avg_wT_err_30=avg_wT_err_30,\
    #invRo_30=invRo_30, avg_uh_30=avg_uh_30, vol_frac_30=vol_frac_30,\
    #eInvRo_100=eInvRo_100, avg_wT_100=avg_wT_100,\
    #eInvRo_err_100=eInvRo_err_100, avg_wT_err_100=avg_wT_err_100,\
    #invRo_100=invRo_100, avg_uh_100=avg_uh_100, vol_frac_100=vol_frac_100)

Nz = 192
Ny = 768
Nx = 768

# loading data from npz file created on the cluster
fn = "vort_frac.npz"
npzfile = np.load(fn)
eInvRo_30=npzfile["eInvRo_30"]
avg_wT_30=npzfile["avg_wT_30"]
eInvRo_err_30=npzfile["eInvRo_err_30"]
avg_wT_err_30=npzfile["avg_wT_err_30"]
invRo_30=npzfile["invRo_30"]
horiz_urms_30=npzfile["horiz_urms_30"]/np.sqrt(Nz*Ny*Nx)
vol_frac_30=npzfile["vol_frac_30"]
eInvRo_100=npzfile["eInvRo_100"]
avg_wT_100=npzfile["avg_wT_100"]
eInvRo_err_100=npzfile["eInvRo_err_100"]
avg_wT_err_100=npzfile["avg_wT_err_100"]
invRo_100=npzfile["invRo_100"]
horiz_urms_100=npzfile["horiz_urms_100"]/np.sqrt(Nz*Ny*Nx)
vol_frac_100=npzfile["vol_frac_100"]

wt_in_turb_100=npzfile["wt_in_turb_100"]
wt_in_turb_30=npzfile["wt_in_turb_30"]
urms_in_turb_100=npzfile["urms_in_turb_100"]
urms_in_turb_30=npzfile["urms_in_turb_30"]


wtMax = [np.max(avg_wT_30), np.max(avg_wT_100)]

fig, ax = plt.subplots(1, 3,figsize=(16,9))
# rearrange ax arary to be indexed 1:4
#ax = [axitem for axitem in [axset for axset in ax]] 
dc = cvar/num_files[0]
col = np.linspace(0+dc,1-dc,len(invRo_30), True)

def ccol(c, dc, n):
    return np.linspace(c-dc,c+dc,n, True)

artists = []

# plots wT B = 30 computation
for i in range(num_files[0]):
    idx_list = np.where(eInvRo_30[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[i],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[0].errorbar(eInvRo_30[i,idx], avg_wT_30[i,idx],
        avg_wT_err_30[i,idx], eInvRo_err_30[i,idx],\
        'none', cm.nipy_spectral_r(clist[j]))
        if j==int(ni/2):
            artists.append(item)

ax[0].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[0].tick_params(axis='both', labelsize=fs-dfs)
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_title(r'$Fr = 0.18$', **font)
fl0 = ax[0].legend(artists,labels_30,title=r'$wT$',loc='lower left',fontsize=fs-dfs)

artists = []

# plots wT B100 computation
for i in range(num_files[1]):
    idx_list = np.where(eInvRo_100[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[cor_idx[i]],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[1].errorbar(eInvRo_100[i,idx], avg_wT_100[i,idx],
        avg_wT_err_100[i,idx], eInvRo_err_100[i,idx],\
        'none', cm.nipy_spectral_r(clist[j]))
        #'none', cm.nipy_spectral_r(colors[cor_idx[i]][idx]))
        if j==int(ni/2):
            artists.append(item)

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
ax[1].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[1].tick_params(axis='both', labelsize=fs-dfs)
ax[1].sharey(ax[0])
ax[1].sharex(ax[0])
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_title(r'$Fr = 0.1$', **font)
fl1 = ax[1].legend(artists,labels_100,title=r'$wT$',loc='lower left',fontsize=fs-dfs)

artists = []

# B = 30 volume fraction data
for i in range(sim_num_files[0]):
    idx_list = np.where(horiz_urms_30[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[i],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[0].scatter(invRo_30[i]/horiz_urms_30[i,idx],
        vol_frac_30[i,idx], s=ms,
        color=cm.nipy_spectral_r(clist[j]))
        #color=cm.nipy_spectral_r(colors[sim_30_idx[i]][2])))
        if j == int(len(idx_list)/2):
            artists.append(item) 

ax[0].legend(artists,labels_30,title=r'$S$',loc='upper left',fontsize=fs-dfs)

artists = []

# B = 100 volume fraction data
for i in range(sim_num_files[1]):
    idx_list = np.where(horiz_urms_100[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[cor_idx[i]],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[1].scatter(invRo_100[i]/horiz_urms_100[i,idx],
        vol_frac_100[i,idx], s=ms,
        color=cm.nipy_spectral_r(clist[j]))
        #color=cm.nipy_spectral_r(colors[sim_30_idx[i]][2])))
        if j == int(ni/2):
            artists.append(item) 

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
#ax[1].set_xlabel(r'$Ro_e^{-1}$', **font)
#ax[1].tick_params(axis='both', labelsize=fs-dfs)
#ax[1].sharey(ax[1,0])
#ax[1].sharex(ax[0,0])
#ax[1].set_yscale('log')
#ax[1].set_xscale('log')
#ax[1].set_title(r'$Fr = 0.1$', **font)
ax[1].legend(artists,labels_100,title=r'$S$',loc='upper left',fontsize=fs-dfs)

artists = []
# wt_rms v u_rms 
for i in range(sim_num_files[0]):
    idx_list = np.where(horiz_urms_30[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[i],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[2].scatter(invRo_30[i]/urms_in_turb_30[i,idx],
        wt_in_turb_30[i,idx],s = ms,c = cm.nipy_spectral_r(clist[j]),
        marker='o')
        if j==int(ni/2):
            artists.append(item)

fl2 = ax[2].legend(artists,labels_30,title=r'$Fr = 0.18$',loc='lower left',fontsize=fs-dfs)

artists = []

for i in range(sim_num_files[1]):
    idx_list = np.where(horiz_urms_100[i,:] > 0)
    ni = len(idx_list[0])
    clist = ccol(col[cor_idx[i]],dc,ni)
    for j, idx in enumerate(idx_list[0]):
        item = ax[2].scatter(invRo_100[i]/urms_in_turb_100[i,idx],
        wt_in_turb_100[i,idx],s = ms,c = cm.nipy_spectral_r(clist[j]),
        marker='^')
        if j==int(ni/2):
            artists.append(item)


ax[2].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[2].tick_params(axis='both', labelsize=fs-dfs)
ax[2].set_yscale('log')
ax[2].set_xscale('log')
ax[2].sharey(ax[0])
ax[2].set_title(r'$wT_{turb}$')
ax[2].legend(artists,labels_100,title=r'$Fr = 0.1$',loc='lower right',fontsize=fs-dfs)


ax[0].add_artist(fl0)
ax[1].add_artist(fl1)
ax[2].add_artist(fl2)
#for i in range(sim_num_files[0]):
    #artists.append(ax[1,0].scatter(invRo_30[i]/avg_uh_30[i],
    #1-vol_frac_30[i], s=ms,
    #color=cm.nipy_spectral_r(colors[sim_30_idx[i]][2])))

#ax[1,0].set_xlabel(r'$Ro_e^{-1}$', **font)
#ax[1,0].set_ylabel(r'$S$',rotation=0, **font)
#ax[1,0].tick_params(axis='both', labelsize=fs-dfs)
#ax[1,0].sharex(ax[0,0])
#ax[1,0].set_yscale('log')
#ax[1,0].set_xscale('log')
#ax[1,0].set_title(r'$Fr = 0.18$', **font)
#ax[1,0].legend(artists,labels_30,loc='lower left',fontsize=fs-dfs)

#artists = []

#for i in range(sim_num_files[1]):
    #artists.append(ax[1,1].scatter(invRo_100[i]/avg_uh_100[i],
    #1-vol_frac_100[i], s=ms,
    #color=cm.nipy_spectral_r(colors[cor_idx[i]][2])))

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
#ax[1,1].set_xlabel(r'$Ro_e^{-1}$', **font)
#ax[1,1].tick_params(axis='both', labelsize=fs-dfs)
#ax[1,1].sharey(ax[1,0])
#ax[1,1].sharex(ax[0,0])
#ax[1,1].set_yscale('log')
#ax[1,1].set_xscale('log')
#ax[1,1].set_title(r'$Fr = 0.1$', **font)
#ax[1,1].legend(artists,labels_100,loc='lower left',fontsize=fs-dfs)

#fig.tight_layout()
fig.subplots_adjust(hspace=.4)
#plt.show()
plt.savefig('flux_vol_frac_plot.pdf')


