
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
from scipy.stats import lognorm,norm
import pandas as pd

#updating the rates of events
num_rates=9; num_states=4; r=np.zeros(num_rates); T=np.zeros([num_rates,num_states])

#function that keeps track of the rate and transition matrices
def update_rates(X,t):
    
    S=X[0]; L=X[1]; A=X[2]; V=X[3]
    
    p=np.random.poisson(pi)
    
    r[0] = aS;             T[0][:]=[1,0,0,0];  #constant production 
    r[1] = dS*S;           T[1][:]=[-1,0,0,0]  #density dependent susceptible death
    r[2] = tau*Bt*S*V;     T[2][:]=[-1,1,0,-1] #latent infection
    r[3] = (1-tau)*Bt*S*V; T[3][:]=[-1,0,1,-1] #active infection
    r[4] = dA*A;           T[4][:]=[0,0,-1,p]  #infected cell burst
    r[5] = gam*V;          T[5][:]=[0,0,0,-1]  #density dependent viral clearance
    r[6] = xi*L;           T[6][:]=[0,-1,1,0]  #latent to active
    r[7] = dL*L;           T[7][:]=[0,-1,0,0]  #latent death
    r[8] = aL*L;           T[8][:]=[0,1,0,0]   #homeostatic division proliferation
    
    
    return r,T #updating the rates of events


#function that solves stochastically using tau-leap method
def simulate_tauleap(t,X0):

    dt=t[1]; x=X0; y=[] #initialize
    for ti in t:
        
        y.append(x) #the list of states

        r,T = update_rates(x,ti) #make new rate vector
        
        E = np.random.poisson(r*dt) #calculate events
        
        dx = np.sum(np.matrix.transpose(T)*E,1) #calculate change in state
        
        x=x+dx #update state variable
        
        x[x<1]=0 #make sure no negative numbers or fractions
        
    return np.array(y)


#parameters for simulations
#cellular parameters
thL = 5.2e-4      # net clearance rate of latent cells
aL  = 0.015;      # proliferation rate of latent cells
dA  = 1.0         # death rate of activated cells
xi  = 0.0001;     # activation rate from latency
aS  = 100         # constant growth rate of suseceptible cells
dS  = 0.03        # death rate of suseceptible cells
Bt  = 1e-4*0.05   # infection rate of T-cells including probability of productively infected
dI  = 1           # death rate of infected cells
pi  = 1e4         # burst rate of virus from cells
gam = 23          # virus clearance rate
mu  = 0.01        # mutation rate
tau = 1e-4        # probability of latency given infction

dL=aL-thL-xi #death rate
l=1-(1+xi/thL)*tau #latency factor

#equilibrium solutions
Seq=gam*dA/Bt0/pi/l
Leq=tau/thL*(gam*dS*dA/Bt0/pi/l-aS)
Aeq=aS*l/dA-gam*dS/Bt0/pi
Veq=aS*pi*l/gam/dA-dS/Bt0
Xeq=np.array([Seq,Leq,Aeq,Veq])

#basic reproductive number
R0=aS*Bt0*pi/gam/dS/dA
print('R_0',R0)
print('ep_c',1-1/R0)

t=np.linspace(0,3*7,1e4) #time in days

#solve the model
tlp_sol = simulate_tauleap(t,X0=np.array([aS/dS,0,0,0]))

#plot the solution
plt.figure(figsize=(6,4),dpi=600)
plt.semilogy(t/7,tlp_sol/5e3)
plt.ylabel('Viral load (copies/mL)')
plt.semilogy(t/7,np.ones(len(t))*Veq/5e3,ls='--',lw=2,color='k')
plt.xlabel('time (weeks)')
plt.tight_layout()

