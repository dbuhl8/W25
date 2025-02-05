import numpy as np
import matplotlib.pyplot as plt


files = ['<cat Om0.5B30Re600Pe60/OUT*',\
    '<cat Om1B30Re600Pe60/OUT*','<cat Om2B30Re600Pe60/OUT*', \
    '<cat Om3B30Re600Pe60/OUT*','<cat Om4B30Re600Pe60/OUT*',\
    '<cat Om5B30Re600Pe60/OUT*','<catOm10B30Re600Pe60/OUT*']
invRo = [0.5, 1, 2, 3, 4, 5, 10]

num_files = len(files)
num_samples = 10

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
        if not line.startswith('#'): #skips headers
            #save data from this line (can choose which data to save) 
            line_elements = [float(x) for x in line.split()]
            np.append(t,line_elements[1])
            np.append(urms,line_elements[3])
            np.append(wT,line_elements[9])
    # perform averaging over eddie turnover times (roughly)
    dt = 2*(np.pi**2)/urms[0]
    idx = np.array([0])
    for j in range(10):
        idx = np.where(t[idx[-1]] <= t <= t[idx[-1]]+dt)
        eInvRo[i,j] = np.sum(urms[idx])*invRo[i]/len(idx)
        eInvRo_err[i,j] = np.std(urms[idx])
        avg_wT[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
    i += 1

fig, ax = plt.subplots(7)

for i in range(num_files):
    ax[i].errorbar(eInvRo[i,:], avg_wT[i,:], avg_wT_err[i,:], eInvRo_err[i,:]) #finish this later
