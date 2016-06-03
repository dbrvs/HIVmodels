#!/usr/bin/python
# Filename: sanctuary_inflamed.py
# ode simulations of the HIV model on ART (R0<1) with reservoir heterogeneity and a sanctuary that has R0>1 and decays in size with time
# updated by DBR 6/2016 #

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

################################################################################
#odes with the tracker    
Yc   = np.zeros(14);             #vector for ODE solutions 
def phaseODEtracker(Xc,t,tau,xi,th,bt,n,a,d,pie,Is_0,S0,phi_s,zeta):
            
    I1    = Xc[0];                              #productively infected cells
    I2    = Xc[1];                              #preintegration cells
    I3em  = Xc[2];  I3cm  = Xc[3]; I3n = Xc[4]; #reservoir cells of different phenotypes
    
    #tracking variables, i for infection, p for proliferation
    I1i   = Xc[5]; 
    I2i   = Xc[6];  I2p   = Xc[7];
    I3iem = Xc[8];  I3pem = Xc[9];
    I3icm = Xc[10]; I3pcm = Xc[11];
    I3in  = Xc[12]; I3pn  = Xc[13];
    
    S=S0*np.exp(-zeta*t);          #susceptible cells
    V=n*(I1+phi_s*Is_0*np.exp(-zeta*t))  #total virus
    
    #set of odes with reservoir heterogeneity
    Yc[0] = tau[0]*bt*S*V + th[0]*I1 + xi[0]*I2 + xi[1]*(I3em+I3cm+I3n)  #I1
    Yc[1] = tau[1]*bt*S*V + th[1]*I2                                     #I2
    Yc[2] = pie[0]*tau[2]*bt*S*V + th[2]*I3em                            #I3em
    Yc[3] = pie[1]*tau[2]*bt*S*V + th[2]*I3cm                            #I3cm
    Yc[4] = pie[2]*tau[2]*bt*S*V + th[2]*I3n                             #I3n

    #equations for tracking variables 
    Yc[5] = tau[0]*bt*S*V + th[0]*I1i                    #I1 infected
    Yc[6] = tau[1]*bt*S*V - d[0]*I2i - xi[0]*I2i         #I2 infected
    Yc[7] = a[0]*I2 - d[0]*I2p - xi[0]*I2p               #I2 proliferating
    Yc[8] = pie[0]*tau[2]*bt*S*V - (d[1] + xi[1])*I3iem  #I3em infected
    Yc[10]= pie[1]*tau[2]*bt*S*V - (d[2] + xi[1])*I3icm  #I3cm infected
    Yc[12]= pie[2]*tau[2]*bt*S*V - (d[3] + xi[1])*I3in   #I3n infected
    Yc[9] = a[1]*I3em - (d[1]+xi[1])*I3pem               #I3em proliferating    
    Yc[11]= a[2]*I3cm - (d[2]+xi[1])*I3pcm               #I3cm proliferating    
    Yc[13]= a[3]*I3n  - (d[3]+xi[1])*I3pn                #I3n proliferating
              
    return Yc
    

################################################################################
#function that evaluates the model
def evaluate_model(t,I1_0,I2_0,I3_0,Is_0,S0,
                    tau,bt,bt_s,n,th,am,xi,d,pie,phi_s,zeta):
    
    Xc0 = [ I1_0,
            I2_0,
            I3_0*pie[0],
            I3_0*pie[1],
            I3_0*pie[2],
            I1_0,
            I2_0,           #I2_0 from activation
            0,              #I2_0 from proliferation
            I3_0*pie[0],
            0,
            I3_0*pie[1],
            0,
            I3_0*pie[2],
            0]
    
    Is_t = Is_0*phi_s*np.exp(-zeta*t) #inf cells in sanctuary over time  
    
    x=odeint(phaseODEtracker, Xc0,t, (tau,xi,th,bt,n,am,d,pie,Is_0,S0,phi_s,zeta) )
    
    #observed infected cells in all compartments from tracking equations
    obs_I=x[:,5]+x[:,6]+x[:,8]+x[:,10]+x[:,12];     
    obs_P=x[:,7]+x[:,9]+x[:,11]+x[:,13]

    #de novo infected cells in all compartments (factors of dt cancel in ratio)
    dno_I=bt*S0*np.exp(-zeta*t)*n*x[:,0]; 
    dno_P=x[:,1]*am[0]+x[:,2]*am[1]+x[:,3]*am[2]+x[:,4]*am[3]
    
    #de novo infected cells in the reservoir I3
    obs3_I=tau[2]*bt*S0*np.exp(-zeta*t)*n*x[:,0]; 
    obs3_P=x[:,2]*am[1]+x[:,3]*am[2]+x[:,4]*am[3]
    
    pct_obs  = (obs_I+Is_t)/(obs_I+obs_P+Is_t)*100                  #observed all
    pct_dno  = (dno_I+Is_t)/(dno_I+dno_P+Is_t)*100                  #de novo all
    pct3_obs = (obs3_I+tau[2]*Is_t)/(obs3_I+obs3_P+tau[2]*Is_t)*100 #de novo I3
    
    I1_t=x[:,0]
    I2_t=x[:,1]
    I3_t=x[:,2]+x[:,3]+x[:,4]
    
    V_t =1e3*n*(x[:,0]+Is_t)  
    
    return I1_t, I2_t, I3_t, Is_t, V_t, pct_obs, pct_dno, pct3_obs               
################################################################################
#plotting function for virus from sanctuary and rest in blood
def main():
        
    tau  = [1,0.1,3e-5]              #the y-intercepts from sarah palmer
    bt0  = 2e-4                      #natural infectivity of HIV [uL/virions-day]
    p    = 1e3                       #[virion/cell-day]
    g    = 23                        #[1/day]
    aS   = 150                       #birth rate susceptible [cells/uL-day]
    dS   = 0.2                       #death rate of susceptible [1/day]
    eps  = 0.95                      #ART efficacy
    s    = [0.46, 0.05, 5e-4]        #multiphasic decay constants [1/day]
    am   = [0.047,0.047,0.015,0.002] #proliferation rates of I2, I3em, I3cm, I3n 
    n    = p/g                       #virions per cell []
    bt   = bt0*(1-eps)               #infectivity on ART

    #global initial conditions
    S0   = 750                       #initial susceptible [cells/uL]
    V0   = 1e5                       #initial virus when ART starts [copies/mL]
    I1_0 = V0/n/1e3;                 #initial infected I1 [cells/uL]
    I2_0 = I1_0/10;                  #initial infected I2 [cells/uL]
    I3_0 = 0.2;                      #initial infected I3 [cells/uL]

    #sanctuary parameters
    eps_s   = 0                      #ART efficacy in the sanctuary
    bt_s    = bt0*(1-eps_s)          #infectivity in the sanctuary
    decay_s = 7.0/365                #decrease in activation pct [1/day] 
    pct_sanc= 1e-5

    #parameters fit to data (decays, initial sanc size, activation rates, death rates)
    th   = [-tau[0]*bt*S0*n-s[0],-s[1],-s[2]]                   

    Is_0 = aS/-th[0]-dS/bt_s/n #initial size of sanctuary 1 in 10^5 

    xi   = [0.08,3e-4]

    d    = [am[0]-th[1]-xi[0],
            am[1]-th[2]-xi[1],
            am[2]-th[2]-xi[1],
            am[3]-th[2]-xi[1]] 

    pie=[0.5,0.4,0.1] #percent em,cm,n respectively
    
    t   = np.linspace(0,365,1e3) #1 year series in days
    
    sol = evaluate_model(t,I1_0,I2_0,I3_0,Is_0,S0,
                    tau,bt,bt_s,n,th,am,xi,d,pie,pct_sanc,decay_s)
                    
    plt.figure(figsize=(8,6))
    
    plt.subplot(211)
    plt.semilogy(t/7, sol[0]*n*1e3, '--',color='dodgerblue', lw=3) #V1
    plt.semilogy(t/7, sol[3]*n*1e3, '--',color='lightcoral', lw=3) #Vs
    plt.semilogy(t/7, sol[4], '-',color='purple', lw=3,alpha=0.5)  #V total
    plt.xlim([0,50])
    plt.xlabel('Time (weeks)')
    plt.ylabel('Viral load (copies/mL)')
    plt.legend([r'$V_1$',r'$V_s$',r'$V=V_1+V_s$'],fontsize=10)
        
    plt.subplot(212)
    plt.plot(t/7,sol[5])
    plt.xlim([0,50])
    plt.xlabel('Time (weeks)')
    plt.ylabel('Observed % infection')
        
    plt.gcf().set_tight_layout(True)
    plt.show()

################################################################################
### MAIN CODE ###
if __name__ == "__main__":
    main()
