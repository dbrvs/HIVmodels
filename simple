#!/usr/bin/python
# 2 simple models of HIV dynamics (inspired by Perelson et al 1996). These models use ordinary differential equations to track the populations of susceptible (S), actively infected and virus producing (A), and virus themselves (V)
# The parameters are taken from the literature, and are checked at the end using the basic reproductive number (the average number of infected cells that arise from a single infected cell in a totally susceptible population) which should be around 8-11 for HIV
# There is a regular susceptible cell growth model, as well as a logistic growth model
# We also check the adiabatic approximation of using the actively infected cells to model virus. The viral timescales are so much faster that it suffices to multiply the infected cells by a constant (the burst size). See the plot for visual confirmation.
# by DBR 7/2015 #

import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

# little function that plots the data (use funny colors!)
def plotSAV(t,S,A,V):
    plt.figure(figsize=(12,6))
    plt.subplot(121)
    plt.semilogy(t, V*1e3, color='maroon', linewidth=3.0,alpha=0.5)
    plt.semilogy(t, A*1e3*pi/gam, color='green', linewidth=1.0,linestyle='--')
    plt.ylabel('Viral concentration (copies per mL)')
    plt.xlabel('Time (days post infection)')
    plt.legend(['ODE solution','calculate V=nA'],loc=4)
    
    plt.subplot(122)
    plt.semilogy(t, S, color='violet', linewidth=3.0)
    plt.semilogy(t, A, color='green', linewidth=3.0)
    plt.legend(['Susceptible','Infected'],loc=4)
    plt.ylabel(r'Concentration (copies per $\mu$L)')
    plt.xlabel('Time (days post infection)')

    plt.show()
    
# parameters for ODE model
aS  = 350    # constant growth rate of suseceptible cells
dS  = 0.3    # death rate of suseceptible cells
Bt  = 1e-4   # infection rate of T-cells
dA  = 1      # death rate of infected cells
pi  = 1e3    # burst rate of virus from cells
gam = 18     # virus clearance rate
T_x = 1e6    # max number of t-cells for logistic growth

tPts = 1e5   #total number of timepoints for ODE solver
tF   = 30    #days

X0   = (1e3,1,0) #initial conditions [S0  A0  V0] cells per uL
    

#set of diff eqs that describes immunologic dynamics
Y=np.zeros((3)); #intialize ode solution vector globally
def SAVmodel(X,t):          
    Y[0] = aS - dS*X[0] - Bt*X[0]*X[2];
    Y[1] = Bt*X[0]*X[2] - dA*X[1];
    Y[2] = pi*X[1] - gam*X[2];
    return Y   # for odeint

def SAVmodel_logistic(X,t):  
    a_t  = aS*(1-(X[0]+X[1])/T_x);  #logistic growth term
    Y[0] = a_t - dS*X[0] - Bt*X[0]*X[2];
    Y[1] = Bt*X[0]*X[2] - dA*X[1];
    Y[2] = pi*X[1] - gam*X[2];
    return Y   # for odeint

tt  = np.linspace(0,tF,tPts);

#RES = spi.odeint(SAVmodel,X0,tt) #typical SAV dynamics

RES = spi.odeint(SAVmodel_logistic,X0,tt) #if want to solve logistic instead

R0=aS*pi*Bt/(dA*dS*gam) #calculate basic reproducive number (~8-11 for HIV)

print 'Basic reproductive number is', R0
   
plotSAV(tt,RES[:,0],RES[:,1],RES[:,2]) #call the plotter function
    

