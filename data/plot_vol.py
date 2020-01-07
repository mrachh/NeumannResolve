from numpy import *
from pylab import *

nlatv = 300
x = loadtxt('helm_neu_vol_'+str(nlatv)+'.dat')
y = loadtxt('plot12.dat2')
xylim = loadtxt('xylim.dat')


xe = abs(x[:,0]-x[:,1])
xe[xe<10**-16] = 10**-16
xe = log10(xe)
xe[x[:,2]<0.9] = inf 
xe = xe.reshape(nlatv,nlatv).transpose()
figure(1)
clf()
imshow(xe,extent=xylim)
cb = colorbar(shrink=0.9,ticks=[-16,-12,-8,-4,0])
cb.ax.set_yticklabels(['$10^{-16}$','$10^{-12}$','$10^{-8}$','$10^{-4}$','$10^{0}$'])
plot(y[:,0],y[:,1],'k-',linewidth=1)
axis('off')
savefig('neu_vol_helm.pdf',bbox_inches='tight')
