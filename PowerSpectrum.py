import Eul_point
import numpy as np
import matplotlib.pyplot as plt



#now: create power spectrum 

delta = 0.001 # time step
c = 0.01          # noise level
ntrial = 10    # nr of spectra to average over (the spectrum is a random variable)
nt = 20000
spec = np.zeros((nt//2+1))  # pre-allocate the spectrum array (frequencies 0 up to nt/2)


for i in range(ntrial): # loop over trials
    print("Starting trial %d..." % (i))
    out = Eul_point.EulerMaruyama(Eul_point.x0, Eul_point.par,nt,delta,c)  # time-step
    dum = np.reshape(out,(nt))                       # for some reason, fft input should be of shape (nt) not (nt,1)
    dum = np.fft.fft(dum)                                  # compute fft
    spec += np.abs(dum[0:nt//2+1])**2       # add for averaging
spec /= ntrial                                                   # divide by nr of trials
freq = np.arange(nt//2+1) / (nt * delta * (Eul_point.t_scale/1000.))  # compute the frequency in Hz from the discrete frequency
plt.loglog(freq,spec)
plt.show()

File = open("spectrum","a")
for i in range(nt//2+1):
    File.write(("%f  %f") % (freq[i],spec[i]))
    File.write("\n")





























