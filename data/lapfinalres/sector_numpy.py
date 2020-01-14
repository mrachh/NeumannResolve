#!/usr/bin/env python

import ccm_widgets as cw
import numpy as np
from matplotlib import cm
cw.init_electron()

nlev = 40
nlat = 300
#fname = 'targ_'+str(nlat)+'.dat'
#rtdat = np.loadtxt(fname)

#rr = np.log10(rtdat[:,1].reshape(nlat,nlat))+14 
#tt= rtdat[:,2].reshape(nlat,nlat)
fname1 = 'scat_'+str(nlat)+'_'+str(nlev)+'.dat'
pref = np.loadtxt(fname1)
ttmin = -0.4776162260822098
ttmax = 0.34374594733735681

zz1 = np.abs(pref[:,0].reshape(nlat,nlat)-pref[:,2].reshape(nlat,nlat))
zz1[zz1<10**-16] = 10**-16
zz1 = np.log10(zz1)
zz1 = np.flipud(zz1)

print(np.min(zz1[:]))
print(np.max(zz1[:]))

r0 = 0.068659

y = cm.viridis(range(256))
X = cw.SectorPlot(
    data_samples=np.array(zz1),
    theta_range=[ttmin, ttmax],
    data_range=[-16,0],
    colorbarticks=[-16,-12,-8,-4,0],
    colorbarticklabels=["10^-16","10^-12","10^-8","10^-4","10^0"],
    colormap=y.T,
    rticks=[r0,0.25+r0,0.5+r0,0.75+r0],
    rticklabels=["10^-12","10^-9","10^-6","10^-3"],
    fontSize=18
)
X.show()
