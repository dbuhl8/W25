import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import torch 

# defining ode functions
def lorenz63(y):
    # y is a numpy array, t is a scalar value
    sigma = 10
    beta = 8./3
    rho = 28
    return np.array([sigma*(y[1]-y[0]),\
        y[0]*(rho-y[2]) - y[1],\
        y[0]*y[1] - beta*y[2]])

def noisylorenz63(y, var):
    noise = np.random.normal(0,var,3)
    return lorenz63(y) + noise


# defining rk4 integrator
def rk4(func, y, dt):
    k1 = func(y)
    k2 = func(y + dt*k1/2)
    k3 = func(y + dt*k2/2)
    k4 = func(y + dt*k3)
    return y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)


def fd1(num,dx):
    # centered finite differnce for intermediate points
    mat = np.zeros([num,num])
    mat[0:num-1,1:num] += np.diag(np.diag(np.ones([num-1,num-1])))
    mat[1:num,0:num-1] -= np.diag(np.diag(np.ones([num-1,num-1])))
    # forward / backward finite difference for start/end points
    mat[0,0] = -2
    mat[0,1] = 2
    mat[num-1,num-2] = -2
    mat[num-1,num-1] = 2
    return mat/(2*dx)


def fd2(num,dx):
    # centered finite difference for intermediate points
    mat = -2*np.diag(np.ones([num,num]))
    mat[0:num-1,1:num] += np.diag(np.ones([num-1,num-1]))
    mat[1:num,0:num-1] += np.diag(np.ones([num-1,num-1]))
    # forward backward difference for start / end points
    mat[0,0] = 1
    mat[0,1] = -2
    mat[0,2] = 1
    mat[num-1,num-3] = 1
    mat[num-1,num-2] = -2
    mat[num-1,num-1] = 1
    return mat/(dx**2)


def gen_feature_matrix(y, num):
    mat = np.zeros([num, 10])
    for i in range(num):
        mat[i, :] = np.array([1, y[i,0], y[i,1], y[i,2], y[i,0]**2, y[i,0]*y[i,1],
                    y[i,0]*y[i,2], y[i,1]**2, y[i,1]*y[i,2], y[i,2]**2])
    return mat

num_points = 10001

t = np.linspace(0, 50, num_points, True)
y = np.zeros([num_points, 3])
y0 = np.array([1, 1, 1])
error = torch.from_numpy(np.zeros_like(y))
y[0, :] = y0
tstep = t[1]-t[0]

# generate time series
for i in range(num_points-1):
    y[i+1, :] = rk4(lorenz63, y[i, :], tstep)

yt = torch.from_numpy(y)
noise = np.array([0.01, 0.1, 0.5])
num_iter = 1
tol = 1e-4

print()
print('-----------------------------------------------')
print()


for i in range(3):
    # generate noisy derivative data
    dy = torch.from_numpy(np.array([noisylorenz63(y[j,:], noise[i])
        for j in range(num_points)]))
    # perform SINDy algorithm
    # define feature matrix and feature coefficients
    phi = torch.from_numpy(gen_feature_matrix(y,num_points))
    theta = torch.linalg.lstsq(phi, dy)[0]
    for j in range(num_iter):
        for k in range(3):
            big_indices = torch.abs(theta[:,k]) >= tol
            small_indices = torch.abs(theta[:,k]) < tol
            theta[small_indices,k] = 0
            theta[big_indices,k] = torch.linalg.lstsq(phi[:,big_indices],dy[:,k])[0]

    print(theta)
    print()
    print('-----------------------------------------------')
    print()
    # perform finite difference method and compare to learned derivatives
    fd = torch.from_numpy(fd1(num_points,tstep))
    thing = torch.matmul(fd,yt) - torch.matmul(phi, theta)
    error[:,i] = torch.sum(torch.abs(thing),1)

fig, ax = plt.subplots()
#colors = cm.rainbow(np.linspace(.5,.7,3))
colors = cm.rainbow([.1, .3, .9])

label = ['Var = 0.01', 'Var = 0.1', 'Var = 0.5']

for i in range(3):
    # plot the error for each solution
    ax.plot(error[:,i], color=colors[i],label=label[i])

ax.set_yscale('log')
ax.legend()
plt.show()
