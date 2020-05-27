# Simple but inefficient way to determine a good number of bins for the frequency plot of intervals.
# In: list z of intervals, minimal number Nmin of data in each bin.
# No output, just produces the frequency plot on a linear and logarithmic scale (the latter to detect power-law decay).
import matplotlib.pyplot as plt
import numpy as np

def makeHistogram(z,Nmin):
    # Determine the minimal and maximal values in the data set and the number of data:
    min = np.amin(z)
    max = np.amax(z)
    Ndat = np.shape(z)[0]
    # Start with 1 bin containing all data:
    Nbins = 1
    least = Ndat
    # Increase the number of bins in a loop, determine the minimal data count, stop when is drops below Nmin:
    while least >= Nmin:
        count = np.zeros((Nbins))            # initialize the bins
        bins = np.linspace(min,max,Nbins+1)  # set the boundaries of the bins (once we have enough data, a multiplicative scaling may be better)
        delt = (max - min) / float(Nbins)    # set the bin width
        # Loop over data, determine the bin they are in and update the "count".
        for i in range(Ndat):
            ind = int(np.amin([np.floor( (z[i] - min) / delt ),Nbins-1]))  # note that the "amin" is necessary to catch the one data point equal to "max"
            count[ind] += 1
        least = np.amin(count)               # find the smallest number of data over all bins
        Nbins += 1
        print("Ndat=%d Nbin=%d least=%d" % (Ndat,Nbins-1,least))
    Nbins = np.amax([Nbins-2,1])             # this is last last value of Nbins for which the minimal data count was greater than or equal to Nmin
    # Re-compute the count with this number of bins (I wrote at the top this is not very efficient!)
    count = np.zeros((Nbins))
    bins = np.linspace(min,max,Nbins+1)
    delt = (max - min) / float(Nbins)
    for i in range(Ndat):
        ind = int(np.amin([np.floor( (z[i] - min) / delt ),Nbins-1]))
        count[ind] += 1
    count = count / (np.sum(count) * delt)   # normalize as a probability density function
    # Plot both on a linear and a logarithmic scale:
    xs = (bins[0:Nbins] + bins[1:Nbins+1]) / 2.0
    plt.plot(xs,count,'-*')
    plt.show()
    plt.loglog(xs,count,'-*')
    plt.show()
    
    
