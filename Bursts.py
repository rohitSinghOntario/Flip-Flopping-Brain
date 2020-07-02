
"""

Simple burst model

periodic spikes of 1 are randomly interupted by periods of 0s.
0s last for an integer # of time periods, based on algebraic prob dist.

"""

import numpy as np
import matplotlib.pyplot as plt

#first defining probability distribution
#letting m run from 1 to 10

#prob dist = rho
rho = []

#power for algebraic prob dist. taken to be 2.75
beta = 2.75

#computing normalization constant
s = 0 
for m in range(1,11):
    s += 1/m**beta
A = 1/s #normalization constant

#creating prob dist.
for m in range(1,11):
    rho.append(A/m**beta)#must be normalized


'''
Function generating periodic bursts, 
skipping integer m # of bursts,
m being an integer drawn from an algebraic distribution

IN: step (timestep), delta (duration of burst), Del (separation of bursts), rho (probability ), T (total duration)

OUT: X

'''
def generateBursts(step, delta, Del, rho, T):
    #time starts at 0
    t = 0
    
    X = []
    
    while t<T:
        #adding peak
        count = 0#counter
        while count<delta:
            X.append(1)
            count+=step
            
        #adding zero
        m = np.random.choice(np.arange(1, 11), p=rho)#random integer from 1 - 10, based on algebraic distribution
        count = 0#counter
        while count<m*(delta + Del):
            X.append(0)
            count+=step
        
        #updating time
        t += delta + m*(delta + Del)
    
    return X




###
#testing
X = generateBursts(0.001,1,10,rho,1000)
plt.title('Bursts')
plt.plot(np.linspace(0,1000,len(X)), X)
plt.show()
###

#now compute power spectrum







    






















