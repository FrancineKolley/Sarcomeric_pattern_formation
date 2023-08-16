% Agent-based simulation of minimal model I
% Copyright (C) <2023>  <Francine Kolley>
% actin scaffold: infitiley long parallel actin lanes
% corresponding to our mean-field simulation assuming steady states = 1
% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% latest version 07_2023

clear all
close all
%% Model parameters --------------------------------------------------------
beta = 1; % [1/s] 
mu=0.85; % [] coupling factor for binding rate for M
zeta=0.85; %  [] coupling factor for binding rate for Z
alpha=1; % []
D=0.01;% Diffusion [m^2/s]
gamma=1; %

%% Numerical parameters ----------------------------------------------------
% numerical grid
dx = 0.02; % [m] bin size
L_sys = 10; % [m] system size
xlist = 0:dx:L_sys-dx; % [m] left boundaries of bins
nx = length(xlist); % [] number of spatial bins

dt = 0.01; % [s] time step
tmax =150; % [s] maximal simulation time
Mast =1; % [1/m] steady-state  concentration myosin
Zast =1; % 1/m] steady-state  concentration myosin


%% translation values
tau=1; %[1/s]
lambda=10000/L_sys; % number density parameter
number_to_concentration = 1/(lambda*dx); % [mol/m]
concentration_to_number = 1/number_to_concentration; % [m/mol]
unit_concentration = 1; % [mol/m]

% number of experiments 
n_exp_max=10;


%% Diffusion probabilty
fun=@(x) 1/sqrt(4*pi*D*dt)*exp(-x.^2/(4*D*dt));
probability=integral(fun,dx/2,Inf);


%% Loop over several experiments
for n_exp=1:n_exp_max

    %% initial condition -------------------------------------------------------
    Mhat=[]; % clear values
    Zhat=[]; % clear vlues
    Mhat = poissrnd( Mast*concentration_to_number, 1, nx ); 
    Zhat = poissrnd( Zast*concentration_to_number, 1, nx );

    %% loop over time-steps
    b=1; %counting time 
    for t=0:dt:tmax

        %% unbinding -----------------------------------------------------------
        n_M = round( sum( Mhat ) ); % total number of bound molecules
        n_Z = round( sum( Zhat ) ); % total number of bound molecules
        n_unbind0M = beta*dt* n_M; % expected number of molecules that unbind in this time step
        n_unbind0Z = beta*dt* n_Z; % expected number of molecules that unbind in this time step
        n_unbindM = poissrnd( n_unbind0M ); % actual number of molecules that unbind in this time step
        n_unbindZ = poissrnd( n_unbind0Z ); % actual number of molecules that unbind in this time step

        n_unbind_M = min(n_unbindM, n_M); % ensure that bin count is always non-negative
        n_unbind_Z = min(n_unbindZ, n_Z); % ensure that bin count is always non-negative
        ind_unbindM = randperm( n_M, n_unbind_M); % indices of molecules that unbind
        ind_unbindZ= randperm( n_Z, n_unbind_Z); % indices of molecules that unbind

        % unbinding process 
        for i = 1:n_unbind_M
            ibin = find( ind_unbindM(i) <= cumsum( Mhat ), 1);
            Mhat(ibin) = Mhat(ibin) - 1;
            if Mhat(ibin) < 0 
                warning('M<0: this should not have happended.\n')
            end
        end
        for i = 1:n_unbind_Z
            ibin = find( ind_unbindZ(i) <= cumsum( Zhat ), 1);
            Zhat(ibin) =Zhat(ibin) - 1;
            if Zhat(ibin) < 0 
                warning('M<0: this should not have happended.\n')
            end
        end
        
        %% binding process -------------------------------------------------
        M_T = Kernel_Titin(1, Mhat, nx,dx)*number_to_concentration;
        
        mast = unit_concentration; % [mol/m]
        zast = unit_concentration; % [mol/m]
        m = number_to_concentration * Kernel_v1(Mhat,nx,dx); % [mol/m] ; note:  ... = M_eps / (lambda*dx)
        z = number_to_concentration * Kernel_v1(Zhat,nx,dx); % [mol/m] ; note:  ... = M_eps / (lambda*dx)
        
        % binding term from mean-field model
        dmdt_bind = beta * mast * exp( mu* ( ( m/mast ) - ( z/zast) - gamma* ( (m/mast) - 1 ).^2 ) ); % [mol/(m s)]
        dzdt_bind = beta * zast * exp(zeta*(alpha*(M_T/mast)+(z/zast)-1-alpha-gamma*((z/zast)-1).^2));
     
        % expected number of molecules that bind in each bin in this time step
        n_bind0M = (concentration_to_number * dmdt_bind * dt); % []
        n_bind0Z = (concentration_to_number * dzdt_bind *dt) ; %[] 
      
        % actual number of molecules that bind in each bin in this time step   
        for l=1:length(n_bind0M)
        n_bind_M = poissrnd( n_bind0M(l) );
        n_bind_Z = poissrnd( n_bind0Z(l) );
        Mhat(l) = Mhat(l) + n_bind_M;
        Zhat(l) = Zhat(l) + n_bind_Z;
        end

      %% Diffusion
     % #1 calculate how many Molecules will jump to left/right bin
     jump_M=floor(Mhat*probability);
     jump_Z=floor(Zhat*probability);
     % #2 substract the jumping molecules from the actual bin
     Mhat=Mhat-2*jump_M;
     Zhat=Zhat-2*jump_Z;
     % #3 update the additonal new molecules for each bin, periodic b.c.
     for k=1:nx
         if k==1
            Mhat(k)=Mhat(k)+jump_M(k+1)+jump_M(nx);
            Zhat(k)=Zhat(k)+jump_Z(k+1)+jump_Z(nx);
         elseif k==nx
             Mhat(k)=Mhat(k)+jump_M(k-1)+jump_M(1);
             Zhat(k)=Zhat(k)+jump_Z(k-1)+jump_Z(1);
         else
             Mhat(k)=Mhat(k)+jump_M(k+1)+jump_M(k-1);
             Zhat(k)=Zhat(k)+jump_Z(k+1)+jump_Z(k-1);
         end
     end 
     
%save for selected time-points
if t<=20 && mod(t,2)<10e-10
   Mhat_save(n_exp,b,:)=Mhat;
   Zhat_save(n_exp,b,:)=Zhat;
   b=b+1;
elseif t>20 && mod(t,10)<10e-10
    Mhat_save(n_exp,b,:)=Mhat;
    Zhat_save(n_exp,b,:)=Zhat;
    b=b+1;
else 
end 


end
end %nexp

