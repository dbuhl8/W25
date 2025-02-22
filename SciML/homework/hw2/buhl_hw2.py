import numpy as np
import matplotlib.pyplot as plt
import torch
import matplotlib as mpl
import matplotlib.cm as cm
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor

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

# need to scramble the data and then allocate portions of the trajectory as
# training and testing data
# scramble
y = np.random.shuffle(y)
num_train = 8000
num_test = 2001

# train model on this trajectory (need to use torch for this)
device = torch.accelerator.current_accelerator().type if torch.accelerator.is_available() else "cpu"
print(f"Using {device} device")

test_data = torch.from_numpy(y[num_train:,:]).to(device)
training_data = torch.from_numpy(y[:num_train,:]).to(device)
batch_size = 64

class lM(torch.nn.Module):
    def __init__(self):
        super(lM, self).__init__()
        self.linear1 = torch.nn.Linear(100, 200)
        self.activation = torch.nn.ReLU()
        self.linear2 = torch.nn.Linear(200, 10)
        self.softmax = torch.nn.Softmax()

    def forward(self, x):
        x = self.linear1(x)
        x = self.activation(x)
        x = self.linear2(x)
        x = self.softmax(x)
        return x

M = lM()

# Create data loaders.
train_dataloader = DataLoader(training_data, batch_size=batch_size)
test_dataloader = DataLoader(test_data, batch_size=batch_size)

loss_fn = nn.MSELoss()
optimizer = torch.optim.SGD(M.parameters(), lr=1e-3)

def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (X, y) in enumerate(dataloader):
        X, y = X.to(device), y.to(device)

        # Compute prediction error
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()

        if batch % 100 == 0:
            loss, current = loss.item(), (batch + 1) * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

def test(dataloader, model, loss_fn):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    model.eval()
    test_loss, correct = 0, 0
    with torch.no_grad():
        for X, y in dataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y).type(torch.float).sum().item()
    test_loss /= num_batches
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")

epochs = 5
for t in range(epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_dataloader, model, loss_fn, optimizer)
    test(test_dataloader, model, loss_fn)
print("Done!")

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
