import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl 

plt.rcParams['text.usetex'] = True

fs = 22
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

invRo_30 = [0.5, 1, 2, 3, 4, 5, 8, 10]
invRo_100 = [0.5, 1, 2, 3, 4, 8]

cor_idx = [0, 1, 2, 3, 4, 6]

labels_100 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=8$']
labels_30 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=5$',r'$\Omega=8$',r'$\Omega=10$']

num_files = [len(files_30), len(files_100)]
num_samples = 5 # num of samples from each timeseries
cvar = 0.2

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
    file.close()

fig, ax = plt.subplots(1, 2,sharex=True,sharey=True,figsize=(16,9))
#fig.suptitle(r'Temp. Flux v $Ro_e^{-1}$', **font)
dc = cvar/num_files[0]
col = np.linspace(0+dc,1-dc,num_files[0], True)
colors = [np.linspace(c-dc,c+dc,num_samples, True) for c in col]

artists = []

for i in range(num_files[0]):
    for j in range(num_samples):
        item = ax[0].errorbar(eInvRo_30[i,j], avg_wT_30[i,j],
        avg_wT_err_30[i,j], eInvRo_err_30[i,j],\
        'none', cm.nipy_spectral_r(colors[i][j]))
        if j==2:
            artists.append(item)

ax[0].set_xlabel(r'$Ro_e^{-1}$', **font)
ax[0].set_ylabel(r'$wT$',rotation=0, **font)
ax[0].tick_params(axis='both', labelsize=fs-5)
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_title(r'$Fr = 0.18$', **font)
ax[0].legend(artists,labels_30,loc='lower left',fontsize=22)


#dc = cvar/num_files[1]
#col = np.linspace(0+dc,1-dc,num_files[1], True)
#colors = [np.linspace(c-dc,c+dc,num_samples, True) for c in col]
artists = []

for i in range(num_files[1]):
    for j in range(num_samples):
        item = ax[1].errorbar(eInvRo_100[i,j], avg_wT_100[i,j],
        avg_wT_err_100[i,j], eInvRo_err_100[i,j],\
        'none', cm.nipy_spectral_r(colors[cor_idx[i]][j]))
        if j==2:
            artists.append(item)

# if tick labels are wanted on right subplot
#ax[1].tick_params(labelleft=True)
ax[1].set_xlabel(r'$Ro_e^{-1}$', **font)
#ax[1].set_ylabel(r'$wT$',rotation=0,labelpad=30, **font)
#ax[1].set_ylabel(r'$wT$',rotation=0, **font)
ax[1].tick_params(axis='both', labelsize=fs-5)
ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_title(r'$Fr = 0.1$', **font)
ax[1].legend(artists,labels_100,loc='lower left',fontsize=fs)

fig.tight_layout()
#plt.show()
plt.savefig('wt_ro.pdf')
