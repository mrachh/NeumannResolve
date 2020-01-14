from numpy import *
from pylab import *

y = loadtxt('sqwts.dat')
nlev = 80
nlat = 300
fname = 'sigma_scat_'+str(nlat)+'_'+str(nlev)+'.dat'
x = loadtxt(fname)

figure(1)
clf()
loglog(x[:,0],abs(x[:,1]),'k.')
savefig('dens_scat.pdf',bbox_inches='tight')


figure(2)
clf()
e1 = abs(x[:,2]-x[:,1])*y
e2 = abs(x[:,3]-x[:,1])*y

loglog(x[:,0],e1,'k.',label='Corner panel - Gauss Legendre nodes')
loglog(x[:,0],e2,'r.',label='Corner panel - custom dirichlet discretization')
xticks([10**-5,10**-15,10**-30])
ylim((10**(-17),10**(3)))
#legend()
fname = 'dens_scat_err_'+str(nlev)+'.pdf'
savefig(fname,bbox_inches='tight')

