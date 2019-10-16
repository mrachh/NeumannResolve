#!/usr/bin/env python

import ccm_widgets as cw
import numpy as np
cw.init_electron()

nlev = 10
nlat = 300
fname = 'targ_'+str(nlat)+'_ref.dat'
rtdat = np.loadtxt(fname)

rr = np.log10(rtdat[:,1].reshape(nlat,nlat))+14 
tt= rtdat[:,2].reshape(nlat,nlat)
fname1 = 'scat_'+str(nlat)+'_'+str(nlev)+'.dat'
pref = np.loadtxt(fname1)
ttmin = -0.4776162260822098
ttmax = 0.34374594733735681

zz1 = np.abs(pref[:,0].reshape(nlat,nlat)-pref[:,1].reshape(nlat,nlat))
zz1[zz1<10**-16] = 10**-16
zz1 = np.log10(zz1)
zz1 = np.flipud(zz1)

X = cw.SectorPlot(
    data_samples=np.array(zz1),
    theta_range=[ttmin, ttmax],
    data_range="auto"
)
X.show()
