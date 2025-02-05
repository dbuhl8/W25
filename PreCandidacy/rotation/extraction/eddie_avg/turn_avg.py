import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl 


files = ['Om0.5B30Re600Pe60/OUT.dat',\
    'Om1B30Re600Pe60/OUT.dat','Om2B30Re600Pe60/OUT.dat', \
    'Om3B30Re600Pe60/OUT.dat','Om4B30Re600Pe60/OUT.dat',\
    'Om5B30Re600Pe60/OUT.dat','Om10B30Re600Pe60/OUT.dat']
invRo = [0.5, 1, 2, 3, 4, 5, 10]

num_files = len(files)
num_samples = 5

eInvRo = np.zeros([num_files, num_samples])
eInvRo_err = np.zeros_like(eInvRo)
avg_wT = np.zeros_like(eInvRo)
avg_wT_err = np.zeros_like(eInvRo)

i = 0
for fn in files:
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
        eInvRo[i,j] = np.sum(urms[idx])*invRo[i]/len(idx)
        eInvRo_err[i,j] = np.std(urms[idx])
        avg_wT[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
    i += 1
    file.close()

fig, ax = plt.subplots()
colors = cm.rainbow(np.linspace(0,1,num_files*num_samples))

for i in range(num_files):
    for j in range(num_samples):
        ax.errorbar(eInvRo[i,j], avg_wT[i,j], avg_wT_err[i,j], eInvRo_err[i,j],\
        'none', colors[i*num_samples + j])

ax.set_xlabel('1/Ro')
ax.set_ylabel('wT')
ax.set_xscale('log')
fig.tight_layout()
plt.savefig('mixing_einvro.pdf')
