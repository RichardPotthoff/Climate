from math import sin,cos,asin,pi,copysign
sign=lambda x:copysign(1,x)
import numpy as np
from fqs import quartic_roots
#approximation of mollweide theta with 4th degree polynomial. 
#max error 0.6% at lat=45°

#  x0=0
#  x1=1
#  y0=0
#  y1=1
#  dy0=2
#  dy1=0
#  ddy1=0
#  m=np.matrix([[ 1  ,   1*x0    ,     1*x0**2     ,    1*x0**3   ,  1*x0**4  ],
#               [ 1  ,   1*x1    ,     1*x1**2     ,    1*x1**3   ,  1*x1**4  ],
#               [ 0  ,   1       ,     2*x0        ,    3*x0**2   ,  4*x0**3  ],
#               [ 0  ,   1       ,     2*x1        ,    3*x1**2   ,  4*x1**3  ],
#               [ 0  ,   0       ,     2           ,    6*x1      , 12*x1**2  ],
#              ])
#  b=np.matrix([[y0],
#               [y1],
#               [dy0],
#               [dy1],
#               [ddy1],
#              ])
#  a=m**-1*b
#  a=list(reversed(a.transpose().tolist()[0]))

def f_lat(theta):
    theta2=abs(theta*2/pi)
    return asin(copysign((2+((-2+theta2)*theta2)*theta2)*theta2,theta))
    
def f_theta(lat):
    return copysign(quartic_roots([1.0, -2.0, 0.0, 2.0, -sin(abs(lat))])[0,1].real*pi/2,lat)

def lat_lon2x_y(lat,lon=0,R=1):
  theta=f_theta(lat)
  x = 2 * R * lon/pi * cos(theta)
  y =     R * sin(theta)
  return x,y
  
def x_y2lat_lon(x,y,R=1):
  theta=asin(y/(R))
  lat=f_lat(theta)
  lon=x/(2*R*cos(theta)) * pi
  return lat,lon

if __name__=='__main__':
  for lat in (-pi/2+0.000,-pi/4,0,pi/4,pi/2-0.000):
    print(f'lat={lat:0.10f} x,y={lat_lon2x_y(lat)}')
  
