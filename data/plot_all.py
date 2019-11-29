from numpy import *
from pylab import *

x = loadtxt('fort.47')
y = loadtxt('fort.48')

nlat = 300

xx = x[:,0].reshape(nlat,nlat)
yy = y[:,0].reshape(nlat,nlat)

z = abs(xx-yy)
z[z<10**-16] = 10**-16
z[z>10**10] = NaN

z = log10(z)
imshow(z)
colorbar()
clim((-16,2))
