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

def bdf3(num,dx):
    if num < 3:
        print("Error: bdf3 cannot take in less than three data points")
        return null
    else:
        mat = np.zeros([num,num])
        bdf3_factor = 11/(6*dx)
        mat += np.diag(np.ones(num-3)*(-2/11),-3)
        mat += np.diag(np.ones(num-2)*(9/11),-2)
        mat += np.diag(np.ones(num-1)*(-18/11),-1)
        mat += np.diag(np.ones(num))
        mat *= bdf3_factor
        # boundary values (will use forward / backward difference)
        # bdf2 
        mat[2,0:3] = (3/(2*dx))*np.array([1/3, -4/3, 1])
        # centered finite difference
        mat[1,0:3] = (1/(2*dx))*np.array([-1, 0, 1])
        # euler forward
        mat[0,0:2] = (1/dx)*np.array([-1, 1])
        return mat

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

act_theta = torch.from_numpy(np.zeros([10,3]))
act_theta[1,0] = -10
act_theta[2,0] = 10
act_theta[1,1] = 28
act_theta[2,1] = -1
act_theta[6,1] = -1
act_theta[3,2] = -8/3
act_theta[5,2] = 1


npact_theta = np.zeros([10,3])
npact_theta[1,0] = -10
npact_theta[2,0] = 10
npact_theta[1,1] = 28
npact_theta[2,1] = -1
npact_theta[6,1] = -1
npact_theta[3,2] = -8/3
npact_theta[5,2] = 1


num_points = 10001

t = np.linspace(0, 50, num_points, True)
y = np.zeros([num_points, 3])
y0 = np.array([1, 1, 1])
y[0, :] = y0
tstep = t[1]-t[0]

# generate time series
for i in range(num_points-1):
    y[i+1, :] = rk4(lorenz63, y[i, :], tstep)

yt = torch.from_numpy(y)
noise = np.array([0.01, 0.1, 0.5])
tol = 1e-4

num_iter = [1, 2, 4, 8, 16]
iter_size = len(num_iter)
error = np.zeros([3,iter_size])

#print()
#print('-----------------------------------------------')
#print()


for i in range(3):
    # generate noisy derivative data
    dy = torch.from_numpy(np.array([noisylorenz63(y[j,:], noise[i])
        for j in range(num_points)]))
    # perform SINDy algorithm
    # define feature matrix and feature coefficients
    phi = torch.from_numpy(gen_feature_matrix(y,num_points))
    for l in range(iter_size):
        theta = torch.linalg.lstsq(phi, dy)[0]
        for j in range(num_iter[l]):
            for k in range(3):
                big_indices = torch.abs(theta[:,k]) >= tol
                small_indices = torch.abs(theta[:,k]) < tol
                theta[small_indices,k] = 0
                theta[big_indices,k] = torch.linalg.lstsq(phi[:,big_indices],dy[:,k])[0]
        #print(num_iter[l])
        #print(theta)
        error[i,l] = torch.sqrt(torch.sum((act_theta - theta)**2))
        #print()
        #print('-----------------------------------------------')
        #print()
    # perform finite difference method and compare to learned derivatives
    #fd = torch.from_numpy(fd1(num_points,tstep))
    #thing = torch.matmul(fd,yt) - torch.matmul(phi, theta)

#fig, ax = plt.subplots()
colors = cm.rainbow(np.linspace(0,1,iter_size))

label = ['N = 1', ' N = 2', 'N = 4', 'N = 8', 'N = 16']

#for i in range(3):
    # plot the error for each solution
    #ax.plot(error[:,i], color=colors[i],label=label[i])
#for i in range(iter_size):
    #ax.plot(noise,error[:,i],color=colors[i],label=label[i])

#ax.set_yscale('log')
#ax.set_xscale('log') #ax.legend()
#ax.set_xlabel('Variance of Noise')
#ax.set_ylabel('MSE error of learned coefficients')
#ax.set_title('SINDy Error (tol = 1e-4)')
#ax.legend(title='N Iterations')


#fig.tight_layout()
#plt.show()
#plt.savefig('noisy_coeff_err.pdf')


# test out the effect different derivative schemes have on SINDy
mfd1 = fd1(num_points,tstep)
mbdf3 = bdf3(num_points,tstep)
fd1_dy = np.matmul(mfd1,y)
bdf3_dy = np.matmul(mbdf3,y)
# perform SINDy algorithm
# define feature matrix and feature coefficients
phi = gen_feature_matrix(y,num_points)
error = np.zeros([2,iter_size])

for l in range(iter_size):
    theta = np.linalg.lstsq(phi, fd1_dy)[0]
    for j in range(num_iter[l]):
        for k in range(3):
            big_indices = np.abs(theta[:,k]) >= tol
            small_indices = np.abs(theta[:,k]) < tol
            theta[small_indices,k] = 0
            theta[big_indices,k] = np.linalg.lstsq(phi[:,big_indices],fd1_dy[:,k])[0]
    error[0,l] = np.sqrt(np.sum((npact_theta - theta)**2, axis=None))

for l in range(iter_size):
    theta = np.linalg.lstsq(phi, bdf3_dy)[0]
    for j in range(num_iter[l]):
        for k in range(3):
            big_indices = np.abs(theta[:,k]) >= tol
            small_indices = np.abs(theta[:,k]) < tol
            theta[small_indices,k] = 0
            theta[big_indices,k] = np.linalg.lstsq(phi[:,big_indices],bdf3_dy[:,k])[0]
    error[1,l] = np.sqrt(np.sum((npact_theta - theta)**2, axis=None))

fig, ax = plt.subplots()

ax.plot(num_iter, error[0,:], color='red',
linestyle='dashed', marker='o', markersize=8, label='FD1')
ax.plot(num_iter, error[1,:], color='blue',
linestyle='dashed', marker='o', markersize=8, label='BDF3')

ax.set_yscale('log')
ax.set_ylabel('MSE Error of learned coefficients')
ax.legend(title='Numerical Method')
ax.set_title('Comparing Derivative Schemes for SINDy')
ax.set_xscale('log')
ax.set_xlabel('Iterations')

fig.tight_layout()
plt.savefig('numerical_coeff_err.pdf')
