#!/usr/bin/python

# by DBR 7/2015 #
# script to read in and plot the data from the c++ SIR model #

import sys
import numpy as np
import matplotlib.pyplot as plt
import time

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

fn = str(sys.argv[1])+".pdf"

# read the file 
f2 = open('HIVsimulation.txt', 'r')
# read the whole file into a single variable, a list of every row of the file.
lines = f2.readlines()
f2.close()

# initialize some variable to be lists:
t=[]; S=[]; I=[]; V=[];

# scan the rows of the file stored in lines, and put the values into S,I,R variables
for line in lines:
    p = line.split()
    t.append(float(p[0]))
    S.append(float(p[1]))
    I.append(float(p[2]))
    V.append(float(p[3]))

plt.subplot(211)    
S_t,I_t=plt.semilogy(t,S,t,I)
plt.setp(S_t, color=(0.172549019607843, 0.627450980392157, 0.1725490196078431), linewidth=3.0)
plt.setp(I_t, color=(0.839215686274509, 0.152941176470589, 0.1568627450980392), linewidth=3.0)
#plt.xlabel('Time (days post infection)')
plt.ylim([1,max(S)*5])
plt.ylabel('Number of T cells')
plt.legend((S_t, I_t), ('Susceptible', 'Infected'))

plt.subplot(212)    
V_t=plt.semilogy(t,V)
plt.ylim([1,max(V)*5])
plt.setp(V_t, color=(0.580392156862745, 0.403921568627451, 0.7411764705882353), linewidth=3.0)
plt.xlabel('Time (days post infection)')
plt.ylabel('Number of virions')
#plt.legend((S_t, I_t, V_t), ('Susceptible', 'Infected', 'Viruses'))


plt.savefig(fn, bbox_inches='tight')

