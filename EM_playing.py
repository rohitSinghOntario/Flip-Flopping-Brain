"""
Purpose: To determine nature of dwell times (exponential or algebraic)

  
1. compute F(x,paramters) 

2. form euler maruyama iteration
    -play around with c value so that we get a reasonable frequency of transitions
    - detect whether x is in state A or state B with a boolean 
    - store dwell times for state A & B
    - with lots of data (transitions) - we wish to see if the tails are algebraic(fat) or exponential
"""
import scipy
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
pee = -0.46*t_scale

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
x0 = x0_per
 

'''
#function computes F(x,u)
#returns F array with 10 components
'''
def F(x,PAR):
    
    #Euler #
    E = np.exp(1)
    
    #What do S1 and S2 represent?
    S1 = PAR[12]/(1+ np.exp(-np.sqrt(2)*(x[0]-PAR[13])/PAR[14]))
    S2 = PAR[15]/(1+ np.exp(-np.sqrt(2)*(x[1]-PAR[17])/PAR[17]))

    #computing components of F, indexed from 0
    f0 = -1-x[0]+(PAR[0]-x[0])/np.abs(PAR[0]+1)*x[2]+(PAR[1]-x[0])/np.abs(PAR[1]+1)*x[4]
    
    f1 = PAR[3]-x[1]+(PAR[0]-x[1])/np.abs(PAR[0]-PAR[3])*x[6]+(PAR[1]-x[1])/np.abs(PAR[1]-PAR[3])*x[8]
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
    
     
  
# storing successive iterations of x in list, for temporary graphing & testing purposes  
x_list = []
t_list = []

'''
Euler Maruyama Integrator

Takes in: initial state: X, parameters: PAR, 
number of time steps: n_timesteps, step size: delta, noise intensity: c


'''
def EulerMaruyama(X, PAR, n_timesteps, delta, c):
    
    for i in range(n_timesteps):
        #computing F
        f = F(X,PAR)
        #drawing random variable from normal distribution
        w = np.random.normal(0,1)
        #updating x
        X = X + delta*f + c*np.sqrt(delta)*w
        
        #need to detect when a transition is made
        #dwell times we get from (delta)(iterations_since_last_switch)
        
        #------------code below this is all temporary & for testing purposes----------
        
        #lets graph x[0] (or any component) vs. time and see if the output is as expected
        x_list.append(X)
        t_list.append(delta*len(x_list))
        
        



#testing above code by graphing an x component vs time_list
        
#calling integrator
EulerMaruyama(x0, par,10000,0.001,1)

#quick function to get list of all nth components of x, to see how they evolve
def nth_component(lst,n): 
    return [element[n] for element in lst] 

#getting a list of a single x-component
x_component_list = nth_component(x_list,0)
#plotting 
plt.plot(t_list, x_component_list)
plt.xlabel("t")
plt.ylabel("x_comp")
plt.show()
        
        
        
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       