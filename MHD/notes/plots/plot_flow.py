import numpy as np
import matplotlib.pyplot as plt


d = 1
M = 10
C = 1

x = np.linspace(-10, 10, 1000)
z = np.linspace(-d, d, 1000)

XX, ZZ = np.meshgrid(x, z)

u = C*(1 - np.cosh(M*ZZ/d)/np.cosh(M))

fig, ax = plt.subplots(figsize=(16,9))

pc1 = ax.pcolor(XX, ZZ, u, cmap='plasma_r')

fig.colorbar(pc1, ax=ax)
ax.set_xticks([])

plt.savefig('bl_hartmann.pdf')

