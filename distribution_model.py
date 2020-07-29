
"""
1. Modelling p(tau) as sum of gaussian bumps
2. compute characteristic phi(omega)
3.compute spectrum S , do we observe 1/omega^somthing with small omega?
"""

import numpy as np
import matplotlib.pyplot as plt


#defining some constants
a  = 2.9 #alpha
b = 0.8 #beta; play around with this
d = 11.5 #delta; period in stable orbit, without noise
M = 10 #number of peaks - should we change this?
a_0 = 0.4 
s_0 = 0.8

#dealing with probability distribution
#computing normalization constant
s = 0 #intermediary sum 
for k in range(1,M+1):
    s += np.sqrt(2*np.pi)*(a_0/k**a)*(s_0*k**b)
A = 1/s #normalization constant

    
#defining charcteristic function
def charfun(w):
    x = 0
    for k in range(1,M+1):
        x+= (a_0/k**a)*(s_0*k**b)*np.exp(-1j*k*d*w-(w**2)*k**((2*b)/2))
    x*=A#normalization const
    return x

#computing power spectrum
#setting range of frequencies
w = np.linspace(0,50,1000)
S = (1 + charfun(w))/(1 - charfun(w))
S = S.real

plt.loglog(w,S)
plt.loglog(w,1/w**0.7)
for k in range(1,M+1):   
    plt.axvline(x=k*2*np.pi/d)#where peaks should be
plt.show()
    




























