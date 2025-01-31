import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# Field Parameters
B_0 = 1
om = 1
eta = 100
alpha = np.sqrt(0.5*om/eta)
print("Wavelength = ", 2*np.pi/alpha, ", Scale Height = ",
1/alpha)


Nz = 1000
Nt = 200
Lz = 2*np.pi/alpha
Period = 4*np.pi/om

# Creating Domain vars
z = np.linspace(0, Lz, Nz, True)
t = np.linspace(0, Period, Nt, True)
Z, T = np.meshgrid(z,t)

B = np.real(B_0 * np.exp(-alpha*Z)*np.cos(om*T)*(np.cos(alpha*Z) + np.sin(alpha*z)))

# Plotting Solution
#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig, ax = plt.subplots()

#surf = ax.plot_surface(Z, T, B, cmap=mpl.colormaps['seismic'])
surf = ax.pcolor(Z, T, B, cmap=mpl.colormaps['seismic'])

ax.set_xlabel("Z")
ax.set_ylabel("T")
fig.colorbar(surf)
fig.tight_layout()
plt.show()
