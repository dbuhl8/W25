import numpy as np
import matplotlib.pyplot as plt

L = 4*np.pi

x = np.linspace(0, L, 50, True)
y = np.linspace(0, L, 50, True)

X, Y = np.meshgrid(x, y)

U = np.cos(2*np.pi*Y/L)*np.sin(2*np.pi*X/L)
V = -np.cos(2*np.pi*X/L)*np.sin(2*np.pi*Y/L)

fig, ax = plt.subplots()
ax.quiver(X, Y, U, V)
plt.show()

