from numpy import *
from pylab import *

fname = 'equi_eig_neu.dat'

f = open(fname)
nverts = int(f.readline())
nroots = int(f.readline())
n = int(f.readline())
nch = int(f.readline())
a=f.readline().split()
nlat = int(a[0])
nt = int(a[1])

xys = zeros((2,n))
dxys = zeros((2,n))
rhs = zeros((2,n))
soln = zeros((2,n))
qwts = zeros(n)

verts = zeros((2,nverts))
vtmp = zeros((2,nverts+1))

for i in range(nverts):
    a= f.readline().split()
    verts[0,i] = float(a[0])
    verts[1,i] = float(a[1])

vtmp[:,0:nverts] = verts
vtmp[:,nverts] = verts[:,0]

for i in range(n):
    a = f.readline().split()
    xys[0,i] = float(a[0])
    xys[1,i] = float(a[1])
    dxys[0,i] = float(a[2])
    dxys[1,i] = float(a[3])
    qwts[i] = float(a[4])

figure(1)
clf()
plot(vtmp[0,:],vtmp[1,:],'k-',linewidth=0.5)
plot(xys[0,:],xys[1,:],'g.',markersize=2)
axis('equal')
show()

arclength = zeros(n)
for i in range(n-1):
    arclength[i+1] = arclength[i] + qwts[i]

ixys = zeros(nch)
ks = zeros(nch)
iscorn = zeros(nch)


for i in range(nch):
    a = f.readline().split()
    ixys[i] = int(a[0])
    ks[i] = int(a[1])
    iscorn[i] = int(a[2])

a = f.readline().split()
xmin = float(a[0])
ymin = float(a[1])
xmax = float(a[2])
ymax = float(a[3])

pot = zeros((2,nt))
targ = zeros((2,nt))
rrt = zeros(nt)
isin = zeros(nt)
errs = zeros(nt)

for irr in range(nroots):
    zk = float(f.readline())
    print(zk)

    rr1 = 0 
    rr2 = 0
    for i in range(n):
        a = f.readline().split()
        rhs[0,i] = float(a[0])
        rhs[1,i] = float(a[1])

    for i in range(nt):
        a = f.readline().split()
        x = float(a[0])
        y = float(a[1])
        targ[0,i] = x 
        targ[1,i] = y
        rrt[i] = sqrt(targ[0,i]**2 + targ[1,i]**2)
        pot[0,i] = float(a[2])
        pot[1,i] = float(a[3])
        isin[i] = int(a[4])
        if(isin[i]>0):
            rr1 = rr1 + pot[0,i]**2
    rr1 = sqrt(rr1)

    for i in range(nt):
        pot[0,i] = pot[0,i]/rr1 
        

    figure(4+irr)
    clf()
    xx = pot[0,:]
    xx[isin<0] = NaN
    xx = xx.reshape(nlat,nlat)
    imshow(xx,extent=[xmin,xmax,ymin,ymax])
    plot(vtmp[0,:],vtmp[1,:],'w-',linewidth=0.5)
    colorbar()
    show()
 
 
