from numpy import *
from pylab import *

nlat = 300

for i in range(1):
    nlev = (i+1)*10
    fname = 'targ_'+str(nlat)+'_ref.dat'
    rtdat = loadtxt(fname)

    rr = log10(rtdat[:,1].reshape(nlat,nlat))+14 
    tt= rtdat[:,2].reshape(nlat,nlat)
    fname1 = 'scat_'+str(nlat)+'_'+str(nlev)+'.dat'
    pref = loadtxt(fname1)

    figure(1) 
    cla()
    clf()
    fig, ax = subplots(subplot_kw=dict(projection='polar'))
    zz1 = abs(pref[:,0].reshape(nlat,nlat)-pref[:,1].reshape(nlat,nlat))
    zz1[zz1<10**-16] = 10**-16
    zz1 = log10(zz1)
    ll2 = linspace(-16,0,100)
    ll = [-16,-12,-8,-4,0]
    llab = ['$10^{-16}$','$10^{-12}$','$10^{-8}$','$10^{-4}$','$10^{0}$']



    cs = ax.contourf(tt, rr, zz1, levels=ll2, vmin=-16.0,vmax=0)
    cbar = fig.colorbar(cs)
    cbar.set_ticks(ll)
    cbar.set_ticklabels(llab)


    ttmin = -0.47761622608220988*360.0/2/pi
    ttmax = 0.34374594733735681*360.0/2/pi

    ax.set_thetamin(ttmin)
    ax.set_thetamax(ttmax)
    ax.set_rlabel_position(ttmax)
    ax.set_rticks([1,4,7,10,13])
    ax.set_yticklabels(['$10^{-13}$','$10^{-10}$','$10^{-7}$','$10^{-4}$','$10^{-1}$'])

    fname2 = 'scat_spdis_err_'+str(nlat)+'_'+str(nlev)+'.pdf'
    savefig(fname2,bbox_inches='tight')

    zz2 = zz1




    figure(1) 
    cla()
    clf()
    fig, ax = subplots(subplot_kw=dict(projection='polar'))
    zz1 = abs(pref[:,0].reshape(nlat,nlat)-pref[:,2].reshape(nlat,nlat))
    zz1[zz1<10**-16] = 10**-16
    zz1 = log10(zz1)
    cs = ax.contourf(tt, rr, zz1, levels=ll2, vmin=-16.0,vmax=0)
    cbar = fig.colorbar(cs)
    cbar.set_ticks(ll)
    cbar.set_ticklabels(llab)
    ttmin = -0.47761622608220988*360.0/2/pi
    ttmax = 0.34374594733735681*360.0/2/pi
    ax.set_thetamin(ttmin)
    ax.set_thetamax(ttmax)
    ax.set_rlabel_position(ttmax)
    ax.set_rticks([1,4,7,10,13])
    ax.set_yticklabels(['$10^{-13}$','$10^{-10}$','$10^{-7}$','$10^{-4}$','$10^{-1}$'])
    fname2 = 'scat_gg_err_'+str(nlat)+'_'+str(nlev)+'.pdf'
    savefig(fname2,bbox_inches='tight')

