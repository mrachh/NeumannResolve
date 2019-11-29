from numpy import *
from pylab import *

fname = 'equi_interior_neu.dat'

f = open(fname)
nverts = int(f.readline())
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
    rhs[0,i] = float(a[5])
    rhs[1,i] = float(a[6])
    soln[0,i] = float(a[7])
    soln[1,i] = float(a[8])

figure(1)
clf()
plot(vtmp[0,:],vtmp[1,:],'k-',linewidth=0.5)
plot(xys[0,:],xys[1,:],'g.',markersize=2)
axis('equal')
show()

arclength = zeros(n)
for i in range(n-1):
    arclength[i+1] = arclength[i] + qwts[i]

figure(2)
clf()
plot(arclength,rhs[0,:],'k.',label='rhs real part',markersize=5)
plot(arclength,rhs[1,:],'kx',label='rhs imag part',markersize=5)
legend()

figure(3)
semilogy(arclength,abs(soln[0,:]),'b.',label='soln real part',markersize=5)
semilogy(arclength,abs(soln[1,:]),'bx',label='soln imag part',markersize=5)
legend()
show()

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

for i in range(nt):
    a = f.readline().split()
    targ[0,i] = float(a[0])
    targ[1,i] = float(a[1])
    pot[0,i] = float(a[2])
    pot[1,i] = float(a[3])

figure(4)
clf()

xx = pot[0,:].reshape(nlat,nlat)
imshow(xx,extent=[xmin,xmax,ymin,ymax])
plot(vtmp[0,:],vtmp[1,:],'w-',linewidth=0.5)
colorbar()
show()
 
y = loadtxt('fort.48')
yy = y[:,0].reshape(nlat,nlat)
figure(5)
clf()

z = abs(xx-yy)
z[z<10**-16] = 10**-16
z = log10(z)
imshow(z,extent=[xmin,xmax,ymin,ymax])
plot(vtmp[0,:],vtmp[1,:],'w-',linewidth=0.5)
colorbar()
show()
