"""
Euler Maruyama iteration, point renewal process
"""

import numpy as np
import matplotlib.pyplot as plt


#Initializing parameters

#excitatory resting potential
her = -70 #mV 
#potential scale
pot_scale = abs(her)

#time scale 
taue = 9 #ms 
t_scale = taue
    
#cell membrane potentials

#excitatory equilibrium potential
hei = 45/pot_scale

#inhibitory equilibrium potential
hii = -90/pot_scale

#tau - inhibitory
taui = 39/t_scale

#inhibitory resting potential
hir = -70/pot_scale

#time constant in the equation for I_ee
at = 0.49*t_scale

#amplitude in RHS of the equation for I_ee
av = 0.81/pot_scale

#time constant in the equation for I_ie
bt = 0.592*t_scale

#amplitude in RHS of the equation for I_ie
bv = 4.85/pot_scale

#external forcing to I_ee (add extra noise here)
#pee = -0.46*t_scale
pee=-4.1

#all other forcings are zero
pei = 0*t_scale
pie = 0*t_scale
pii = 0*t_scale

#product of # of excitatory-excitatory for excitatory-inhibitory connections with amplitude me, nbee*emax
nbee_emax = 3034*0.5*t_scale

#theta_e, threshold value for h_e
te = -50/pot_scale

#steepness of the transition
se = 5/pot_scale

#product of # of inhibitory-excitatory for inhibitory-inhibitory connections with amplitude me, nbie*imax
nbie_imax = 536*0.5*t_scale

#theta_i, threshold value for h_i
ti = -50/pot_scale

#steepness of the transition for h_i
si = 5/pot_scale

#putting all parameters into par
par = np.array([hei,hii,taui,hir,at,av,bt,bv,pee,pei,pie,pii,nbee_emax,te,se,nbie_imax,ti,si])



#defining initial x (got this from corticalEEGmodel)
x0_eq = np.array([-1.0368,-1.0200,0.1344,0.0,0.2,0.0,0.1639,0.0,0.2,0.0])
x0_per = np.array([-0.9598,-0.9305,3.4119,-5.0216,2.9430,-3.5510,3.4414,-5.0216,2.9430,-3.5510])
x0 = x0_eq
 

#precomputing constants for F where possible
sqrt2 = np.sqrt(2)
f0_cons1 = np.abs(par[0]+1)#f0, constant 1
f0_cons2 = np.abs(par[1]+1)

f1_cons1 = np.abs(par[0]-par[3])
f1_cons2 = np.abs(par[1]-par[3])



'''
#function computes F(x,u)
#returns F array with 10 components
'''
def F(x,PAR):
    
    #Euler #
    E = np.exp(1)
    
    #What do S1 and S2 represent?
    S1 = PAR[12]/(1+ np.exp(-sqrt2*(x[0]-PAR[13])/PAR[14]))
    S2 = PAR[15]/(1+ np.exp(-sqrt2*(x[1]-PAR[16])/PAR[17]))

    #computing components of F, indexed from 0
    f0 = -1-x[0]+x[2]*(PAR[0]-x[0])/f0_cons1+x[4]*(PAR[1]-x[0])/f0_cons2
    
    f1 = PAR[3]-x[1]+x[6]*(PAR[0]-x[1])/f1_cons1+x[8]*(PAR[1]-x[1])/f1_cons2
    f1 = f1/PAR[2]
    
    f2 = x[3]
    
    f3 = -2*PAR[4]*x[3]-PAR[4]**2*x[2]+PAR[5]*PAR[4]*E*(S1+PAR[8])
    
    f4 = x[5]
    
    f5 = -2*PAR[6]*x[5]-PAR[6]**2*x[4]+PAR[7]*PAR[6]*E*(S2+PAR[10])
    
    f6 = x[7]
    
    f7 =-2*PAR[4]*x[7]-PAR[4]**2*x[6]+PAR[5]*PAR[4]*E*(S1+PAR[9])
    
    f8 = x[9]
    
    f9 = -2*PAR[6]*x[9]-PAR[6]**2*x[8]+PAR[7]*PAR[6]*E*(S2+PAR[11])
    
    
    return np.array([f0,f1,f2,f3,f4,f5,f6,f7,f8,f9])
    
     
  


#threshold value at which state change is detected (wrt X[0]) 
#threshold was obtained visually - perhaps theres a better way
switch = -0.8


'''
Euler Maruyama Integrator
IN: initial state: X, parameters: PAR, 
number of time steps: n_timesteps, step size: delta, noise intensity: c

OUT:h_e list
'''
def EulerMaruyama(X, PAR, n_timesteps, delta, c):
    #list storing all X
    x_list = []
    
    #to update dwell times before storing
    dwell_time = 0
    
    #precomputing as much as possible
    rand = [np.random.normal(0,1) for i in range(n_timesteps)]
    cSqrtDelta = c*np.sqrt(delta)
    
    #file storing dwell times
    he_file = open("dwell_times.txt", "a")
    
    for i in range(n_timesteps):
        
        #storing 'previous' value of X[0] for comparisons, when detecting a 'switch'
        X_0_prev = X[0]
        
        
        #computing F
        f = F(X,PAR)
        
        #updating x
        X = X + delta*f + cSqrtDelta*rand[i]
        
        #building x list
        x_list.append(X)
        
        
        
    
        #identifying state changes & updating/storing dwell times
        
        #one dwell time = the time between consecutive 'shots up'
        #must identify shots up - then store current dwell time and reset
        
        if X[0]<switch:
            #at base level, so update dwell time
            dwell_time += delta
        elif X_0_prev<switch:
            #this conditional is called if X[0]>=switch && prev X[0]<switch
            #therefore a switch is detected
            
            #store dwell time
            he_file.write(str(dwell_time))
            he_file.write("\n")
            
            #reset dwell time
            dwell_time = 0
            
     
    he_file.close()
    #returning h_e
    return [element[0] for element in x_list] 
