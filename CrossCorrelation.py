# -*- coding: utf-8 -*-
"""
checking correlation between dwell times
"""
import numpy as np
import matplotlib.pyplot as plt

#reading all data from textfile
he_file = open("dwell_times.txt", "r")
stringList = he_file.readlines()
he_file.close()

#converting from strings
dwell_times = [float(x) for x in stringList]


# #of data points
M = len(dwell_times)

#average dwell\
avg = sum(dwell_times)/M

#data indexed for delay
k = range(1,M+1)

#correlation function
c = []

for delay in k:
    ck = 0
    for i in range(1,M-delay):
        ck+=(dwell_times[i]-avg)*(dwell_times[i+delay]-avg)/M
    c.append(ck)

#plotting
plt.title('dwell correlation')
plt.plot(k,c)
plt.show()
        


