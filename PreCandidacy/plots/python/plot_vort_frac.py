import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
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

# there is a section missing from vort_frac.py where data extraction is
# performed 

# this is the command used to save the data
#np.savez("vort_frac", eInvRo_30=eInvRo_30, avg_wT_30=avg_wT_30,\
    #eInvRo_err_30=eInvRo_err_30, avg_wT_err_30=avg_wT_err_30,\
    #invRo_30=invRo_30, avg_uh_30=avg_uh_30, vol_frac_30=vol_frac_30,\
    #eInvRo_100=eInvRo_100, avg_wT_100=avg_wT_100,\
    #eInvRo_err_100=eInvRo_err_100, avg_wT_err_100=avg_wT_err_100,\
    #invRo_100=invRo_100, avg_uh_100=avg_uh_100, vol_frac_100=vol_frac_100)


# loading data from npz file created on the cluster
fn = "vort_frac.npz"
npzfile = np.load(fn)
eInvRo_30=npzfile["eInvRo_30"]
avg_wT_30=npzfile["avg_wT_30"]
eInvRo_err_30=npzfile["eInvRo_err_30"]
avg_wT_err_30=npzfile["avg_wT_err_30"]
invRo_30=npzfile["invRo_30"]
avg_uh_30=npzfile["avg_uh_30"]
vol_frac_30=npzfile["vol_frac_30"]
eInvRo_100=npzfile["eInvRo_100"]
avg_wT_100=npzfile["avg_wT_100"]
eInvRo_err_100=npzfile["eInvRo_err_100"]
avg_wT_err_100=npzfile["avg_wT_err_100"]
invRo_100=npzfile["invRo_100"]
avg_uh_100=npzfile["avg_uh_100"]
vol_frac_100=npzfile["vol_frac_100"]

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

artists = []

for i in range(sim_num_files[0]):
    artists.append(ax[1,0].scatter(invRo_30[i]/avg_uh_30[i],
    vol_frac_30[i], s=ms,
    color=cm.nipy_spectral_r(colors[sim_30_idx[i]][2])))

ax[1,0].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[1,0].set_ylabel(r'$S$',rotation=0, **font)
ax[1,0].tick_params(axis='both', labelsize=fs-dfs)
ax[1,0].set_yscale('log')
ax[1,0].set_xscale('log')
ax[1,0].set_title(r'$Fr = 0.18$', **font)
ax[1,0].legend(artists,labels_30,loc='lower right',fontsize=fs-dfs)

artists = []

for i in range(sim_num_files[1]):
    artists.append(ax[1,1].scatter(invRo_100[i]/avg_uh_100[i],
    vol_frac_100[i], s=ms,
    color=cm.nipy_spectral_r(colors[cor_idx[i]][2])))

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
ax[1,1].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[1,1].tick_params(axis='both', labelsize=fs-dfs)
ax[1,1].sharey(ax[1,0])
ax[1,1].sharex(ax[1,0])
ax[1,1].set_yscale('log')
ax[1,1].set_xscale('log')
ax[1,1].set_title(r'$Fr = 0.1$', **font)
ax[1,1].legend(artists,labels_100,loc='lower right',fontsize=fs-dfs)

#fig.tight_layout()
fig.subplots_adjust(hspace=.4)
#plt.show()
plt.savefig('flux_vol_frac_plot.pdf')


