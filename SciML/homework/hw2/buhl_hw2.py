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

print("\n----------------------------------------------------------------------------\n")
# ------------------------------------------------------------------------------------------
print(" Part a: \n\n")

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
    h, c, b = (.5, 8, 10) # need to check these coefficients
    dxy = np.zeros_like(xy)
    # transcribes the data into smaller arrays
    x = xy[:8]
    y = np.array(
        [xy[8*(i+1):8*(i+1)+8] for i in range(8)])
    # computes the derivative with two loops
    for i in range(8):
        dxy[i] = x[i-1]*(x[(i+1)%8]-x[i-2]) - x[i] + F - (h*c/b)*np.sum(y[i,:])
        for j in range(8):
            dxy[8*(i+1)+j] = -c*b*y[i,(j+1)%8]*(y[i,(j+2)%8] - y[i,j-1]) - c*y[i,j] +\
                            (h*c/b)*x[i]
    return dxy

def rk4(func, y, dt):
    k1 = func(y)
    k2 = func(y + dt*k1/2)
    k3 = func(y + dt*k2/2)
    k4 = func(y + dt*k3)
    return y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

# initializing trajectory (why not a vector of ones)
y0 = np.random.randn(72)
z0 = y0+(1e-5)*np.random.randn(72)
nt = 20000
dt = .005
y = np.zeros([nt+1, 72])
z = np.zeros([nt+1, 72])
y[0,:] = y0
z[0,:] = z0
F = 20


for i in range(nt):
    y[i+1,:] = rk4(lorenz96, y[i,:],dt)
    z[i+1,:] = rk4(lorenz96, z[i,:],dt)

t = np.linspace(0, nt*dt, nt+1,True)

fig1, ax1 = plt.subplots(1, 2)
fig2, ax2 = plt.subplots()

ax1[0].plot(t, np.sqrt(np.sum((y[:,:8] - z[:,:8])**2, axis=1)))
ax1[0].set_yscale("log")
ax1[0].set_xlabel("Time")
ax1[0].set_title("Using a small perturbation")


def sum_y(y):
    ysum = torch.zeros([len(y[:,0]),8]).to(dtype=dtype, device=device)
    for i in range(8):
        ysum[:,i] = torch.sum(y[:,8*(i):8*(i+1)], 1)
    return ysum

# need to scramble the data and then allocate portions of the trajectory as
# training and testing data
# scramble
np.random.shuffle(y)
num_train = 19500
num_test = 501
num_epochs = 10
num_ts = 20

# train model on this trajectory (need to use torch for this)
device = torch.accelerator.current_accelerator().type if torch.accelerator.is_available() else "cpu"
dtype = torch.float32
print(f"Using {device} device")

test_data = torch.from_numpy(y[num_train:,:]).to(dtype=dtype, device=device)
training_data = torch.from_numpy(y[:num_train,:]).to(dtype=dtype, device=device)
batch_size = 64

class lM(torch.nn.Module):
    def __init__(self):
        super(lM, self).__init__()
        self.l1 = torch.nn.Linear(8, 512)
        self.l2 = torch.nn.Linear(512, 512)
        self.l3 = torch.nn.Linear(512, 512)
        self.l4 = torch.nn.Linear(512, 512)
        self.l5 = torch.nn.Linear(512, 512)
        self.l6 = torch.nn.Linear(512, 512)
        self.l7 = torch.nn.Linear(512, 512)
        self.l8 = torch.nn.Linear(512, 512)
        self.l9 = torch.nn.Linear(512, 512)
        self.l10 = torch.nn.Linear(512, 64)
        self.act = torch.nn.ReLU()

    def forward(self, x):
        x = self.l1(x)
        x = self.act(x)
        x = self.l2(x)
        x = self.act(x)
        x = self.l3(x)
        x = self.act(x)
        x = self.l4(x)
        x = self.act(x)
        x = self.l5(x)
        x = self.act(x)
        x = self.l6(x)
        x = self.act(x)
        x = self.l7(x)
        x = self.act(x)
        x = self.l8(x)
        x = self.act(x)
        x = self.l9(x)
        x = self.act(x)
        x = self.l10(x)
        return x

M = lM().to(device)
loss_fn = nn.MSELoss()
optimizer = torch.optim.Adam(M.parameters(), lr=1e-4)

def train(trajectory, model, loss_fn, optimizer):
    global batch_size
    size = len(trajectory[:,0])
    batch_list = [[batch_size*i,batch_size*(i+1)] for i in
    range(int(size/batch_size))]
    batch_list[-1][-1] = size
    #print(batch_list)
    model.train()
    #for batch, (X, y) in enumerate(dataloader):
    for batch, idx  in enumerate(batch_list):
        x = trajectory[idx[0]:idx[1],:8]
        y = trajectory[idx[0]:idx[1],8:]

        # Compute prediction error
        pred = model(x)
        #loss = loss_fn(pred, sum_y(y))
        loss = loss_fn(pred, y) + 100*loss_fn(sum_y(pred), sum_y(y))

        # Backpropagation
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        if batch%20 == 0:
            print("Loss :", loss.item(), ", Batch: [",\
            batch_size*batch,"/",num_train,"]")

def test(trajectory, model, loss_fn):
    size = len(trajectory[:,0])
    global batch_size
    batch_list = [[batch_size*i,batch_size*(i+1)] for i in
    range(int(size/batch_size))]
    batch_list[-1][-1] = size
    #print(batch_list)
    num_batches = len(batch_list)
    #for batch, (X, y) in enumerate(dataloader):

 
    model.eval()
    pred_loss, correct = 0, 0
    sum_loss = 0
    with torch.no_grad():
        for idx  in batch_list:
            x = trajectory[idx[0]:idx[1],:8]
            y = trajectory[idx[0]:idx[1],8:]
            pred = model(x)
            pred_loss += loss_fn(pred, y).item()
            sum_loss += loss_fn(sum_y(pred), sum_y(y)).item()
            #correct += (pred.argmax(1) == y).type(torch.float).sum().item()
    pred_loss /= num_batches
    sum_loss /= num_batches
    print("Avg Y Pred Loss: ", pred_loss,"\n Avg Y Sum Loss: ", sum_loss)

for j in range(num_ts):
    y = np.zeros([nt+1, 72])
    y[0,:] = np.random.randn(72)

    for i in range(nt):
        y[i+1,:] = rk4(lorenz96, y[i, :],dt)
    np.random.shuffle(y)
    test_data = torch.from_numpy(y[num_train:,:]).to(dtype=dtype, device=device)
    training_data = torch.from_numpy(y[:num_train,:]).to(dtype=dtype, device=device)
    for k in range(num_epochs):
        #print(f"Epoch {t+1}\n-------------------------------")
        train(training_data , M, loss_fn, optimizer)
        test(test_data, M, loss_fn)
        print("  ")
        print("Finished Epoch: ", k+1)
        print("  ")
    print(" ------------- ")
    print("  Completed Timeseries: ", j+1)
    print(" ------------- ")

print("Completed Training (F = 20), testing prediction error")

# generating a new test set (trajectory)
y[0,:] = np.random.randn(72)
for i in range(nt):
    y[i+1,:] = rk4(lorenz96, y[i,:],dt)

np.random.shuffle(y)
ytest = torch.from_numpy(y).to(dtype=dtype,device=device)
test_data = torch.from_numpy(y).to(dtype=dtype, device=device)
test(test_data,M,loss_fn)

ts = 100
def plot_sum_np(point):
    return [np.sum(point[8*(i+1):8*(i+2)]) for i in range(8)]
def plot_sum_torch(point):
    return [torch.sum(point[8*(i+1):8*(i+2)]).cpu().detach().numpy() for i in range(8)]
ax2.plot(plot_sum_np(y[ts,:]), 'bo', label=r'Act')
ax2.plot(plot_sum_torch(M(ytest[ts,:8])), 'r+', label=r'Pred')
ax2.set_title(r'Differences between Actual and Prediction')
ax2.set_xlabel('i')
ax2.set_ylabel(r'$\sum_jY_{i,j}$', rotation=0)
ax2.legend()


print("\n----------------------------------------------------------------------------\n")
# ----------------------------------------------------------------------------
print(" Part b: \n\n")

# problem b is to solve the same problem but with F = 24 (not sure what F is)
# psuedo code: (repeat problem a with different parameters)
# 1. generate testing data (with F = 24)
# 2. compute error of predictions and sum value
# 3. retrain the model with new data

# compute new trajectory 
#y0 = np.ones(72)

y = np.zeros([nt+1, 72])
y[0,:] = y0

F = 24
for i in range(nt):
    y[i+1,:] = rk4(lorenz96, y[i, :],dt)

np.random.shuffle(y)
ytest = torch.from_numpy(y).to(dtype=dtype,device=device)
test_data = torch.from_numpy(y[num_train:,:]).to(dtype=dtype, device=device)
training_data = torch.from_numpy(y[:num_train,:]).to(dtype=dtype, device=device)

# compute error in model predictions
print("\n Testing Model on Trajectory with F = 24\n")
test(ytest, M, loss_fn)
print("\n\n")


# retrain the model 
print(" Retraining the Model on a Trajectory with F = 24 ")
epochs = num_epochs
for k in range(epochs):
    train(training_data , M, loss_fn, optimizer)
    test(test_data, M, loss_fn)
    print("  ")
    print("Finished Epoch: ", k+1)
    print("  ")
print("\nCompleted Training (F = 24)\n")


print("\n----------------------------------------------------------------------------\n")
# ------------------------------------------------------------------------------------------
print(" Part c: \n\n")


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
    global dtype
    global device
    global M
    h, c, b = (.5, 8, 10) # need to check these coefficients
    dx = np.zeros_like(x)
    y =  M(torch.from_numpy(x).to(dtype=dtype,device=device))
    # computes the derivative with two loops
    for i in range(8):
        dx[i] = x[i-1]*(x[(i+1)%8]-x[i-2]) - x[i] + F \
                -(h*c/b)*torch.sum(y[8*i:8*(i+1)])
    return dx


x = np.zeros([nt+1,8])
y = np.zeros([nt+1,72])
y[0,:] = np.random.randn(72)
x[0,:] = y[0,:8]

F = 20

for i in range(nt):
    x[i+1,:] = rk4(lorenz96M, x[i,:],dt)
    y[i+1,:] = rk4(lorenz96, y[i,:],dt)

# compute divergence from actual trajectory
div_x = np.sqrt(np.sum((y[1:,:8]-x[1:,:])**2)/nt)
print("RMS Divergence from Actual Trajectory: ", div_x, " (Normalized MSE Loss)")

t = np.linspace(0,nt*dt,nt+1,True)
# visualize the "divergence"

ax1[1].plot(t[1:], np.sqrt(np.sum((y[1:,:8]-x[1:,:])**2, axis=1)))
ax1[1].set_xlabel("Time")
ax1[1].set_yscale('log')
ax1[1].set_title(r'Using the Model v Actual')
ax1[0].set_ylabel("D", rotation=0)
fig1.suptitle(r'Distance between Trajectories')

fig1.tight_layout()

fig1.savefig("traj_div.pdf")
fig2.savefig("y_sum_error.pdf")

torch.save(M.state_dict(), "lorenz96_Y")
# net.load_state_dict(torch.load(load_model_path + 'FNO2D_Eulerstep_MSE_Loss_Randomized_Cahn_Hilliard_modes_12_wavenum_50_lead_10.pt'))

