% ############## Mean-field model #############
% Copyright (C) <2023>  <Francine Kolley>
% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% latest version 07_2023
% numerical simulation for Mean-field model
% euler-scheme - time iteartion
% stencil.m - Laplacian 5point stencil
% Kernel.m - convolution with double gaussian kernel 
% periodic boundary conditions 

clear all 
close all 

%% system parameters 
% parameter which define the size of the grid and space between points
L_sys=10;        % System size
n_x=500;         % number of grid points
N=0:n_x-1;      
N=N*L_sys/n_x;  % whole grid vector
dx=L_sys/n_x;   % Space between grid points 

% Diffusion
Dm=0.01;
Da=0.01;

% non-local-interaction length
l_0=1;
sigma=l_0/5; %width of kernel 

% additional model parameter
mu=0.85;
zeta=0.85;
alpha=1; 
gamma=1;
beta=1;
m_ast=1;
z_ast=1;

%% initial conditions 
conc_M(1,:)=1*ones(1,n_x)+0.01*rand(1,n_x);
conc_Z(1,:)=1*ones(1,n_x)+0.01*rand(1,n_x);


%% numerical iterarion 
% Laplacian (5-point stencil, call function stencil)
[L]=stencil(n_x,dx);

% parameters for time iteration 
tmax=100;
dt=0.01;
for j=0:dt:tmax
  
    % Call function for Convolution
    [M_Z]=Kernel(l_0,conc_M,n_x,dx,sigma); %Titin(l_t,M,n_x,dx);


    conc_M=conc_M+(dt*Dm*(L*(conc_M)'))'+dt*(-beta*conc_M+beta*m_ast*exp(mu*(conc_M-conc_Z-1*gamma*(conc_M-1).^2)));
    conc_Z=conc_Z+(dt*Da*(L*(conc_Z)'))'+dt*(-beta*conc_Z+beta*z_ast*exp(zeta*(conc_Z+alpha*M_Z-m_ast-alpha*z_ast-1*gamma*(conc_Z-1).^2)));


end

figure()
plot(N,conc_M,'b')
hold on
plot(N,conc_Z,'g')
ylim([0 2])
