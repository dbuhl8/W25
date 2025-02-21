import numpy as np
import matplotlib.pyplot as plt
import torch
import matplotlib as mpl
import matplotlib.cm as cm

# the idea behind this assignment is to solve the Lorenz' 96 system

# problem a is to have a model predict all of the Y variables given the X data
# set
# psuedo code:
# 1. generate training data (x,y)
# 2. design a model to compute y from x
# 3. generate testing data (a different trajectory)
# 4. compute error of predictions and sum value

def lorenz96(xy):
    # this takes in xy which is a row vector with the first 8 components being
    # x_k, and the next 64 components y_jk. i.e. k in [1,8], j in [1,64]
    global F 
    h, c, b = (1, 1, 1) # need to check these coefficients
    dxy = np.zeros_like(xy)
    # transcribes the data into smaller arrays
    x = xy[:8]
    y = np.array(
        [xy[8*(i+1):8*(i+1)+8] for i in range(8)])
    # computes the derivative with two loops
    for i in range(8)
        dxy[i] = x[i-1]*(x[i+1]-x[i-2]) - x[i] + F - (h*c/b)*np.sum(y[i,:])
        for j in range(8)
            dxy[8*(i+1)+j] = -c*b*y[i,j+1]*(y[i,j+2] - y[i,j-1]) - c*y[i,j] +\
                            (h*c/b)*x[i]
    return dxy

def rk4(func, y, dt):
    k1 = func(y)
    k2 = func(y + dt*k1/2)
    k3 = func(y + dt*k2/2)
    k4 = func(y + dt*k3)
    return y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

# initializing trajectory (why not a vector of ones)
y0 = np.ones(72)
nt = 10000
dt = .001

y = np.zeros([nt, 72])
y[0,:] = y0

F = 8

for i in range(nt):
    y[i+1,:] = rk4(lorenz96, y[:,i],dt)

# train model on this trajectory (need to use torch for this)

# oop

# problem b is to solve the same problem but with F = 24 (not sure what F is)
# psuedo code: (repeat problem a with different parameters)
# 1. generate testing data (with F = 24)
# 2. compute error of predictions and sum value
# 3. retrain the model with new data

# compute new trajectory 
y0 = np.ones(72)
nt = 10000
dt = .001

y = np.zeros([nt+1, 72])
y[0,:] = y0

F = 24
for i in range(nt):
    y[i+1,:] = rk4(lorenz96, y[:,i],dt)

# compute error in model predictions

# retrain the model 

# problem c is to couple the model to the dynamics of the system, i.e.
# numerically integrate the model where sum_j y_j,k is obtained using the model
# M1.
# pseudo code:
# 1. write an rk4 solver to integrate X using only X and M1(X) to compute 
# Sum_j Y_jk
# 2. Stabilize the system (there might be errors in M1 which propagate into X
# and blow up the system)
# 3. Compute error of the trajectory

def lorenz96M(x):
    # computes the X components of the lorenz 96 system using M to predict y
    global F 
    h, c, b = (1, 1, 1) # need to check these coefficients
    dx = np.zeros_like(x)
    # computes the derivative with two loops
    for i in range(8)
        dx[i] = x[i-1]*(x[i+1]-x[i-2]) - x[i] + F #-(h*c/b)*np.sum(M(x))
    return dx

x0 = np.ones(8)
x = np.zeros([nt+1,72])

F = 24

for i in range(nt):
    x[i+1,:] = rk4(lorenz96M, x[i,:],dt)

# compute divergence from actual trajectory

# attempt to stabilize the system
