from numpy import *
from pylab import *

y = loadtxt('sqwts.dat')

ffs = ['ref','px','rand','scat']
for i in range(4):
    fname = 'dir_sigma_'+ffs[i]+'_100.dat'
    x = loadtxt(fname)
    figure(1)
    clf()
    loglog(x[:,0],abs(x[:,1]),'k.')
    xticks([10**-5,10**-15,10**-30])
    savefig('dens_lap_dir_'+ffs[i]+'.pdf',bbox_inches='tight')
    figure(2)
    clf()
    e = abs(x[:,1]-x[:,2])
    e[e<10**-16] = 10**-16
    e = e*y
    loglog(x[:,0],e,'k.')
    xticks([10**-5,10**-15,10**-30])
    savefig('dens_err_lap_dir_'+ffs[i]+'.pdf',bbox_inches='tight')


nlevs = ['10','20','40','80']
for i in range(4):
    figure(2)
    clf()
    figure(3)
    clf()
    for j in range(4):
        fname = 'sigma_'+ffs[i]+'_300_'+nlevs[j]+'.dat'
        x = loadtxt(fname)

        if(j == 1):
            figure(1)
            clf()
            loglog(x[:,0],abs(x[:,1]),'k.')
            xticks([10**-5,10**-15,10**-30])
            savefig('dens_lap_neu_'+ffs[i]+'.pdf',bbox_inches='tight')

        e = abs(x[:,1]-x[:,3])
        e[e<10**-16] = 10**-16
        e = e*y
        figure(2)
        loglog(x[:,0],e,'.',label=nlevs[j])
        
        figure(4)
        clf()
        loglog(x[:,0],e,'k.',label=nlevs[j])
        xticks([10**-5,10**-15,10**-30])
        savefig('dens_err_lap_neu_'+ffs[i]+'_'+nlevs[j]+'.pdf',bbox_inches='tight')
        

        e = abs(x[:,1]-x[:,2])
        e[e<10**-16] = 10**-16
        e = e*y
        figure(3)
        loglog(x[:,0],e,'.',label=nlevs[j])

        figure(5)
        clf()
        loglog(x[:,0],e,'k.',label=nlevs[j])
        xticks([10**-5,10**-15,10**-30])
        savefig('dens_err_lap_neu_gl_'+ffs[i]+'_'+nlevs[j]+'.pdf',bbox_inches='tight')
        

    figure(2)
    xticks([10**-5,10**-15,10**-30])
    legend()
    savefig('dens_err_lap_neu_'+ffs[i]+'.pdf',bbox_inches='tight')

    figure(3)
    xticks([10**-5,10**-15,10**-30])
    legend()
    savefig('dens_err_lap_neu_gl_'+ffs[i]+'.pdf',bbox_inches='tight')
show()
