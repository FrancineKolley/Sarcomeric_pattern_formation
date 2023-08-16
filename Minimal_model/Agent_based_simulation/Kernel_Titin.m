% Function for the interaction with Titin 
% Copyright (C) <2023>  <Francine Kolley>
% area normalized to 1!
% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% latest version 07_2023

function[M_T]=Kernel_Titin(l_t,A,n_x,dx)
%% variables necessary for A_T
%sigma=l_t/5; %width of kernel 
var=0.04;
xlist_for_At=(-n_x/2+1:1:n_x/2)*dx; %xlist for A_t

%% A_T
%%% Scaling factor 2 for all cases including more than one sarcomere
M_T= (exp(-(xlist_for_At-l_t).^2/(2*var))+ exp(-(xlist_for_At+l_t).^2/(2*var)))/(2*sqrt(2*pi*var));

%% convolution 
xlist_for_concA=(0:n_x-1)*dx;
xlist_for_conv=xlist_for_At(1)+(0:2*n_x-2)*dx;
convolution_At=conv(A,M_T)*dx;

%% periodic boundary conditions 
 tleft = xlist_for_concA(1);
 tright = xlist_for_concA(end);
 jleft = find(xlist_for_conv>=tleft,1);
 jright = find(xlist_for_conv>=tright,1);
 
 xlist_for_convAt_periodic_bc=xlist_for_conv(jleft:jright); %( maybe add +/- 1) looks ok
 M_T = convolution_At(jleft:jright);
 M_T(end-(jleft-2):end) =  M_T(end-(jleft-2):end) + convolution_At(1:jleft-1);
 M_T(1:length(convolution_At)-jright) =  M_T(1:length(convolution_At)-jright) + convolution_At(jright+1:end);
end 