from numpy import *
from pylab import *

ttmin = -0.477616226
ttmax = 0.34374594


x = [0,0.6*cos(ttmin)]
y = [0,0.6*sin(ttmin)]

plot(x,y,'k-')


x = [0,0.6*cos(ttmax)]
y = [0,0.6*sin(ttmax)]

plot(x,y,'k-')

x = [0.6*cos(ttmin),0.9*cos(ttmin)]
y = [0.6*sin(ttmin),0.9*sin(ttmin)]

plot(x,y,'k--')


x = [0.6*cos(ttmax),0.9*cos(ttmax)]
y = [0.6*sin(ttmax),0.9*sin(ttmax)]

plot(x,y,'k--')


x = [0.9*cos(ttmin),cos(ttmin)]
y = [0.9*sin(ttmin),sin(ttmin)]

plot(x,y,'k-')


x = [0.9*cos(ttmax),cos(ttmax)]
y = [0.9*sin(ttmax),sin(ttmax)]


plot(x,y,'k-')



nlev = 4
rticks = 2.0**(-arange(nlev))
print(rticks)
xx = rticks
rn1 = [-sin(ttmin),cos(ttmin)]
rn2 = [-sin(ttmax),cos(ttmax)]
h = 0.005
for i in range(nlev):
    x0 = [xx*cos(ttmin),xx*sin(ttmin)]

    x = [x0[0]-rn1[0]*h,x0[0]+rn1[0]*h]
    y = [x0[1]-rn1[1]*h,x0[1]+rn1[1]*h]
    plot(x,y,'k-')

    x0 = [xx*cos(ttmax),xx*sin(ttmax)]

    x = [x0[0]-rn2[0]*h,x0[0]+rn2[0]*h]
    y = [x0[1]-rn2[1]*h,x0[1]+rn2[1]*h]
    plot(x,y,'k-')


axis('equal')
axis('off')
savefig('grid.pdf',bbox_inches='tight')
