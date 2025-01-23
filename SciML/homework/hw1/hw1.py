import torch
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu"
)
n = 50
xstart = 1
xend = 10
x = np.linspace(xstart,xend,n,True)
y = 2*x + np.cos(25*x) + np.random.normal(0,1,n )

tn = 50
tx = np.sort(np.random.uniform(xstart,xend,tn))
ty = 2*tx + np.cos(25*tx) + np.random.normal(0,1,tn )

iterate = True

if iterate:

    pst = 5
    pen = 500
    pstep = 1

    P = np.arange(pst, pen, pstep)
    err = np.zeros_like(P)
    terr = np.zeros_like(P)

    for i in range(np.size(P)):

        # Building feature matrix: Random Fourier Features
        om = np.random.normal(0,1,P[i]) # number of features is given by p
        phi = np.cos(2*np.pi*np.outer(x, om))
        tphi = np.cos(2*np.pi*np.outer(tx, om))

        # Doing a linear regression (case by case for size of phi)
        if P[i] < n: # do a least squares solution
            TH = np.linalg.lstsq(phi,y)[0]
        else: # do a psuedoinverse solution (can switch to np.pinv)
            TH = np.matmul(np.linalg.pinv(phi),y)
            #U, S, Vh = np.linalg.svd(phi, True, True)
            #TH = np.dot(Vh[0:n,:].conj().T,
            #np.multiply(np.divide(1,S,where=S>0),np.matmul(U.conj().T,y)))

        # compute the error 
        err[i] = np.linalg.norm(y-np.matmul(phi, TH),2)
        terr[i] = np.linalg.norm(ty-np.matmul(tphi, TH),2)

    fig, (ax1, ax2, ax3) = plt.subplots(3,1)
    #fig, ax = plt.subplots()
    #ax.semilogy(P,err, '-', color='orange', linewidth=0.5)
    #ax.semilogy(P, terr, 'b--', linewidth=0.5)
    #ax.title.set_text("Error")
    #ax.set_ylabel(r'$||y - A\Theta^*||_2$')
    #ax.set_xlabel('Number of Features')

    ax1.plot(x,y,'bx',x,np.matmul(phi, TH),'g--',linewidth=1,markersize=3)
    ax1.title.set_text("Training Data")
    ax1.set_ylim([0,20])
    ax2.plot(tx,ty,'bx',tx,np.matmul(tphi, TH),'g--',linewidth=1,markersize=3)
    ax2.title.set_text("Testing Data")
    ax2.set_ylim([0,20])
    ax3.semilogy(P,err, '-',color='orange',linewidth=0.5,label="Train")
    ax3.semilogy(P, terr, 'b--', linewidth=0.5, label="Test")
    ax3.legend(loc="upper right")
    ax3.title.set_text("Error")
    ax3.set_ylabel(r'$||y - A\Theta^*||_2$')
    ax3.set_xlabel('Number of Features')
    #plt.axis([5,100,1E-7,1000])

    fig.tight_layout()

    plt.savefig('hw1s2.pdf') 

else:
    p = 50
    # Building feature matrix: Random Fourier Features
    om = np.random.normal(0,1,p) # number of features is given by p
    phi = np.cos(2*np.pi*np.outer(x, om))
    tphi = np.cos(2*np.pi*np.outer(tx, om))

    # Doing a linear regression (case by case for size of phi)
    if p < n: # do a least squares solution
        TH = np.linalg.lstsq(phi,y)[0]
    else: # do a psuedoinverse solution (can switch to np.pinv)
        TH = np.matmul(np.linalg.pinv(phi,1E-8),y)
        #U, S, Vh = np.linalg.svd(phi, True, True)
        #TH = np.dot(Vh[0:n,:].conj().T,
        #np.multiply(np.divide(1,S,where=S>0),np.matmul(U.conj().T,y)))

    #print(np.divide(1,S,where=S>0))

    # compute the error 
    err = np.linalg.norm(y-np.matmul(phi, TH),2)
    terr = np.linalg.norm(ty-np.matmul(tphi, TH),2)

    fig, (ax1, ax2) = plt.subplots(2,1)

    print("Test error:",err,", Training Error:", terr)

    ax1.plot(x,y,'bx',x,np.matmul(phi, TH),'g--')
    ax2.plot(tx,ty,'bx',tx,np.matmul(tphi, TH),'g--')
    plt.axis([0,10,0,20])
    plt.show()


