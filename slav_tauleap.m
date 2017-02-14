%Matlab script that uses tau-leap to stochastically simulated the SLAV
%model, that's the susceptible, latently infected, actively infected and
%virus model

%% updating the rates of events
function [t,logV_sol]=slav_tauleap()
    
    num_sims=10; %how many curves to make
    
    tF=20; %how many days to sample
    sim_tpts=1e3; %initial number of tpts to simulate (before downsampling below)
    
    t=linspace(0,tF,sim_tpts); %time in days to simulate

    % add Gaussian noise and downsample (uncomment to get no effect)
    noise=0.5; %magnitude of measurement noise in logs?
    %noise=0; %magnitude of measurement noise in logs?
    tpts=100; %number of timepoints after downsampling
    %tpts=sim_tpts; %number of timepoints after downsampling
    
    ds_factor=round(length(t)/tpts); %downsample factor to get right timepoints
    tt=downsample(t,ds_factor); %downsampled time
    simulated_data=zeros(num_sims+1,length(tt)); %matrix for downsampled data
    simulated_data(1,:)=tt; %first column in matrix is timeseries
    
    %simulate a few times
    for i=1:num_sims
        tic; logV_sol=simulate_slav(t); toc;
        vv=downsample(logV_sol+randn(length(logV_sol),1)*noise,ds_factor);
        simulated_data(i+1,:)=vv;
    end
    
figure();
plot(tt,simulated_data(2:end,:),'o')
xlabel('time (days)')
ylabel('log_{10} copies HIV per mL')

save('simulated_data')
    
    
function V_sol = simulate_slav(t)
    %parameters for simulations
    vol = 1e3;           % volume of simulation
    thL = 5.2e-4;      % net clearance rate of latent cells
    aL  = 0.015;       % proliferation rate of latent cells
    dA  = 1.0;         % death rate of activated cells
    xi  = 0.0001;      % activation rate from latency
    aS  = 100*vol;         % constant growth rate of suseceptible cells
    dS  = 0.03;        % death rate of suseceptible cells
    Bt  = 5e-6/vol*10^randn();   % infection rate of T-cells including probability of productively infected
    p   = 1e4;         % burst rate of virus from cells
    gam = 23;          % virus clearance rate
    tau = 1e-4;        % probability of latency given infction
    dL=aL-thL-xi;      % death rate of latent cells

    %equilibrium solutions for dynamical system with virus, can be useful
    %l=1-(1+xi/thL)*tau; %latency factor
    %Seq=gam*dA/Bt/p/l;
    %Leq=tau/thL*(gam*dS*dA/Bt/p/l-aS);
    %Aeq=aS*l/dA-gam*dS/Bt/p;
    %Veq=aS*p*l/gam/dA-dS/Bt;
    %Xeq=[Seq,Leq,Aeq,Veq];
    
    %basic reproductive number
    R0=aS*Bt*p/gam/dS/dA;
    disp(['R0=' num2str(R0)]); %note if R0<1 infection burns out

    %solve the model
    X0=[aS/dS,0,1,0]; %start at disease free equilibrium with a single infected cell
    tlp_sol = simulate_tauleap(t,X0,aS,dS,tau,Bt,dA,p,gam,xi,dL,aL);

    %V_sol=tlp_sol(:,4); %just keep the virus for now
    V_sol=log10(tlp_sol(:,4)+30); %just keep the virus for now, as typical in logs


%%
%function that keeps track of the rate and transition matrices
function [r,T] = update_rates(X,t,r,T,aS,dS,tau,Bt,dA,p,gam,xi,dL,aL)
    
    S=X(1); L=X(2); A=X(3); V=X(4);
    
%%
    ps=poissrnd(p); %stochastic burst size
    r(1) = aS;             T(1,:)=[1,0,0,0];  %constant production 
    r(2) = dS*S;           T(2,:)=[-1,0,0,0];  %density dependent susceptible death
    r(3) = tau*Bt*S*V;     T(3,:)=[-1,1,0,-1]; %latent infection
    r(4) = (1-tau)*Bt*S*V; T(4,:)=[-1,0,1,-1]; %active infection
    r(5) = dA*A;           T(5,:)=[0,0,-1,ps];  %infected cell burst
    r(6) = gam*V;          T(6,:)=[0,0,0,-1];  %density dependent viral clearance
    r(7) = xi*L;           T(7,:)=[0,-1,1,0];  %latent to active
    r(8) = dL*L;           T(8,:)=[0,-1,0,0];  %latent death
    r(9) = aL*L;           T(9,:)=[0,1,0,0];   %homeostatic division proliferation
%%    
    

%function that solves stochastically using tau-leap method
function y=simulate_tauleap(t,X0,aS,dS,tau,Bt,dA,p,gam,xi,dL,aL)

    num_rates=9; num_states=4; r=zeros([num_rates,1]); T=zeros([num_rates,num_states]);

    dt=t(2); x=X0; y=zeros([length(t),num_states]); %initialize
    
    for inx=1:length(t)
                
        y(inx,:)=x; %the list of states

        [r,T]=update_rates(x,t(inx),r,T,aS,dS,tau,Bt,dA,p,gam,xi,dL,aL); %make new rate vector
        
        E = poissrnd(r*dt); %calculate events
        
        dx = T'*E; %calculate change in state
        
        x=x+dx'; %update state variable
        
        x(x<1)=0; %make sure no negative numbers or fractions
        
    end 