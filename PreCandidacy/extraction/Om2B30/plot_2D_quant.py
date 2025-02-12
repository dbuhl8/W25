import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


file = open("Om2B30.dat", 'r')
Nx, Ny, Nt,dx,dy = file.readline().split()
Nx = int(Nx) - 2
Ny = int(Ny) - 2
Nt = int(Nt)
dx = float(dx)
dy = float(dy)

X = np.zeros((Nx, Ny), dtype=np.int32)
Y = np.zeros_like(X)
T = np.zeros(Nt)

Ux = np.zeros((Nx, Ny, Nt), dtype=np.float64)
Uy = np.zeros_like(Ux)
Uz = np.zeros_like(Ux)
Uz_rms = np.zeros_like(Ux)
Wz = np.zeros_like(Ux)
TDisp = np.zeros_like(Ux)
Lz = np.zeros_like(Ux)
MDisp = np.zeros_like(Ux)
Baro_ER = np.zeros_like(Ux)
GR = np.zeros_like(Ux)

# it would be nice if the data_extraction script wrote the number of lines, and
# timesteps contained on the first line of the file
t = 0
k = -1
counter = 0
for line in file:
    line_elements = [float(x) for x in line.split()]
    if len(line_elements) > 0 :
        i = int(line_elements[0] - 1)
        j = int(line_elements[1] - 1)
        X[i, j] = i+1
        Y[i, j] = j+1
        if (t != line_elements[2]):
            t = line_elements[2]
            k = k+1
            T[k] = t
        Ux[i,j,k], Uy[i,j,k], Uz[i,j,k], Wz[i,j,k], TDisp[i,j,k], \
            Uz_rms[i,j,k], Lz[i,j,k], Baro_ER[i,j,k], MDisp[i,j,k], GR[i,j,k] \
            = line_elements[3:]

fig, ((ax1, ax2), (ax3,ax4)) = plt.subplots(2,2)

ax1.pcolor(X,Y,Wz[:,:,-1],cmap=mpl.colormaps['seismic'])
ax2.pcolor(X,Y,Lz[:,:,-1])
ax3.pcolor(X,Y,Uz_rms[:,:,-1])
ax4.pcolor(X,Y,GR[:,:,-1])

fig.tight_layout()
plt.show()
