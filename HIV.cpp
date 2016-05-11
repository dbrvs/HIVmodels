/* ****************************************************************

from K&R code

DBR 2015

******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>

// Set up basic parameters

double alph=0.5;
double beta=1e-7;
double delt=0.5;
double kapp=200;
double gamm=5;
double Tmax=1e6;

double S0=1e6;
double I0=0;
double V0=1;

double tF=70;

// Set up variables and rates of change
double t,S,I,V,Pop[3];
double dPop[3];

// The system of ODEs describing HIV dynamics
void Diff(double Pop[3])
{
  double tmpS, tmpI, tmpV, alphat;
  
  tmpS=Pop[0]; tmpI=Pop[1]; tmpV=Pop[2]; //temp variables for neatness
  
  alphat=alph*(1-(tmpS+tmpI)/Tmax); //logistic growth term
  
  dPop[0] = alphat*tmpS - beta*tmpS*tmpV; // dS/dt
  dPop[1] = beta*tmpS*tmpV - delt*tmpI;   // dI/dt
  dPop[2] = kapp*delt*tmpI - gamm*tmpV;   // dV/dt

  return;
}

// Runge-Kutta 4th order solver
void Runge_Kutta(double step)
{
  int i;
  double dPop1[3], dPop2[3], dPop3[3], dPop4[3];
  double tmpPop[3], initialPop[3];

  initialPop[0]=S; initialPop[1]=I; initialPop[2]=V;

  Diff(initialPop);
  for(i=0;i<3;i++)
    {
      dPop1[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
    }

  Diff(tmpPop);
  for(i=0;i<3;i++)
    {
      dPop2[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop2[i]/2;  
    }

  Diff(tmpPop);
  for(i=0;i<3;i++)
    {
      dPop3[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop3[i]; 
    }

  Diff(tmpPop);

  for(i=0;i<3;i++)
    {
      dPop4[i]=dPop[i];
      tmpPop[i]=initialPop[i]+(dPop1[i]/6+dPop2[i]/3+dPop3[i]/3+dPop4[i]/6)*step;
    }


  S=tmpPop[0]; I=tmpPop[1]; V=tmpPop[2];

  return;
}

// Run the integrator and output the data
int main()
{
  double step=0.01;
  double t=0;
     
  S=S0; I=I0; V=V0; //initial conditions
  
  std::ofstream myfile ("HIVsimulation.txt");
  if (myfile.is_open())
  {
	  //myfile << "----t----S----I----V----" << std::endl; //header line in file
	  
	  do
	  { 
		  myfile << t << " " << S << " " << I << " " << V << " " << std::endl;
	  	  Runge_Kutta(step);  t+=step;
	  } 
	  while(t<tF);
	  
      myfile.close();
   }
   else std::cout << "Unable to open file";
  
  
  return 0;

}