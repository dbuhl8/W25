import numpy as np
import matplotlib.pyplot as plt

ep = 5
x = np.linspace(-ep,ep,100)
y = np.linspace(-ep,ep,100)

XX, YY = np.meshgrid(x,y)
alpha = 1.5
BX = YY
BY = alpha**2*XX

fig, ax = plt.subplots()

ax.streamplot(XX, YY, BX, BY)

fig.tight_layout()
plt.savefig('ex5_fieldlines.pdf')
