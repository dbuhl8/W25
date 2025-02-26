import numpy as np


# these are Finite difference schemes for fields taken from the simdat files
def FDX(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # centered finite difference formula (inside)
    for i in range(nx-2):
        dx_field[:,:,:,i+1] = (0.5/dx)*(field[:,:,:,i+2]-field[:,:,:,i])
    return dx_field

def FDY(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,:,i+1,:] = (0.5/dy)*(field[:,:,i+2,:]-field[:,:,i,:])
    return dy_field

def FDZ(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dz_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dz
    dz_field[:,nz-1,:,:] = (field[:,nz-1,:,:]-field[:,nz-2,:,:])/dz
    # centered finite difference formula (inside)
    for i in range(nz-2):
        dz_field[:,i+1,:,:] = (0.5/dz)*(field[:,i+2,:,:]-field[:,i,:,:])
    return dz_field

def FD4X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,1] = (field[:,:,:,2]-field[:,:,:,0])/(2*dx)
    dx_field[:,:,:,nx-2] = (field[:,:,:,nx-1]-field[:,:,:,nx-3])/(2*dx)
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # 4th order centered finite difference formula (inside)
    for i in range(nx-4):
        dx_field[:,:,:,i+2] = (1/(12*dx))*(-field[:,:,:,i+4]+8*field[:,:,:,i+3]\
                                           -8*field[:,:,:,i+1]+field[:,:,:,i])
    return dx_field

def FD4Y(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # 4th order centered finite difference formula (inside)
    for i in range(ny-4):
        dy_field[:,:,i+2,:] = (1/(12*dy))*(-field[:,:,i+4,:]+8*field[:,:,i+3,:]\
                                           -8*field[:,:,i+1,:]+field[:,:,i,:])
    return dy_field

def FD4Z(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    dz_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dz
    dz_field[:,1,:,:] = (field[:,2,:,:]-field[:,0,:,:])/(2*dz)
    dz_field[:,nz-2,:,:] = (field[:,nz-1,:,:]-field[:,nz-3,:,:])/(2*dz)
    dz_field[:,nz-1,:,:] = (field[:,nz-1,:,:]-field[:,nz-2,:,:])/dz
    # 4th order centered finite difference formula (inside)
    for i in range(nz-4):
        dz_field[:,i+2,:,:] = (1/(12*dz))*(-field[:,i+4,:,:]+8*field[:,i+3,:,:]\
                                           -8*field[:,i+1,:,:]+field[:,i,:,:])
    return dz_field

def FD6X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    dx_field[:,:,:,1] = (field[:,:,:,2]-field[:,:,:,0])/(2*dx)
    dx_field[:,:,:,2] = (-field[:,:,:,4]+8*field[:,:,:,3] \
                        -8*field[:,:,:,1]+field[:,:,:,0])/(12*dx)
    dx_field[:,:,:,nx-3] = (-field[:,:,:,nx-1]+8*field[:,:,:,nx-2]\
                            -8*field[:,:,:,nx-5]+field[:,:,:,nx-6])/(12*dx)
    dx_field[:,:,:,nx-2] = (field[:,:,:,nx-1]-field[:,:,:,nx-3])/(2*dx)
    dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    # 4th order centered finite difference formula (inside)
    for i in range(nx-6):
        dx_field[:,:,:,i+3] = (1/(60*dx))*(-field[:,:,:,i+6]+9*field[:,:,:,i+5]\
                                           -45*field[:,:,:,i+4]\
                                           +45*field[:,:,:,i+2]-9*field[:,:,:,i+1]\
                                           +field[:,:,:,i])
    return -dx_field

def FD6Y(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    dy_field[:,:,2,:] = (-field[:,:,4,:]+8*field[:,:,3,:] \
                        -8*field[:,:,1,:]+field[:,:,0,:])/(12*dy)
    dy_field[:,:,ny-3,:] = (-field[:,:,ny-1,:]+8*field[:,:,ny-2,:]\
                            -8*field[:,:,ny-5,:]+field[:,:,ny-6,:])/(12*dy)
    dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    # 4th order centered finite difference formula (inside)
    for i in range(ny-6):
        dy_field[:,:,i+3,:] = (1/(60*dy))*(-field[:,:,i+6,:]+9*field[:,:,i+5,:]\
                                           -45*field[:,:,i+4,:]\
                                           +45*field[:,:,i+2,:]-9*field[:,:,i+1,:]\
                                           +field[:,:,i,:])
    return -dy_field

