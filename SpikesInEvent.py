
"""
This is the EulPoint code modified to 
count number of spikes per event (set to 15 ms for now)

Intermittent poisson process?

What is the tail behaviour of # of spikes per event?
"""


import numpy as np

#spikes within 'event_length' of eachother are considered one event
event_length = 15 #15ms for now - 1.5 periods

#list storing number of spikes per event (per delta)
spikes = []

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
    
     
  


#threshold value at which spike is detected (wrt X[0]) 

switch = -0.8


'''
Euler Maruyama Integrator
IN: initial state: X, parameters: PAR, 
number of time steps: n_timesteps, step size: delta, noise intensity: c

'''
def EulerMaruyama(X, PAR, n_timesteps, delta, c):
    
    #precomputing as much as possible
    rand = [np.random.normal(0,1) for i in range(n_timesteps)]
    cSqrtDelta = c*np.sqrt(delta)
    rand_adjust = np.array([0,0,0,1,0,0,0,0,0,0])
    
    #file storing spikes per event
    spike_file = open("spikes.txt", "a")
    
    #temp storage of spike # per event
    n = 0
    
    #temp storage of spikes - do not want to write to file every event
    spike_temp = []
    
    #tracking time delay - renews every spike or when reaches zero
    t_delay = event_length
    
    #boolean telling us if we are in an event
    inEvent = False
    
    for i in range(n_timesteps):
        
        #storing 'previous' value of X[0] for comparisons, when detecting a 'switch'
        X_0_prev = X[0]
        
        
        #computing F
        f = F(X,PAR)
        
        #updating x
        X = X + delta*f + cSqrtDelta*rand[i]*rand_adjust
        
        
        #updating t_delay & spike #
        if inEvent == True:
            t_delay -= delta
            if t_delay<=0:
                #reset delay time
                t_delay = event_length
                #temp spike storage
                spike_temp.append(n)  
                #reset n 
                n=0
                #reset bool
                inEvent = False
                
        #periodically writing to file - every 100 events
        if len(spike_temp)>=100:
            #store spike #s in file
            for num in spike_temp:
                spike_file.write(str(num))
                spike_file.write("\n")
            #reset temp spike storage
            spike_temp = []
            
        #identifying spikes/events
        if X_0_prev<switch and X[0]>switch:
            #increase n by one (spike # in one event)
            n+=1
            
            #trigger start of event
            if inEvent == False:
                inEvent = True
            
            #reset delay time
            t_delay = event_length
            
    spike_file.close()
    
     
EulerMaruyama(x0, par, 2000000, 0.001, 0.1)






        

        
        













