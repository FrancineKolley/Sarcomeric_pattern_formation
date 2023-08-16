% Small gaussian kernel to ensure bin coupling
% Copyright (C) <2023>  <Francine Kolley>
function[M_K]=Kernel_v1(M,n_x,dx)
%% variables necessary for A_T
mu=0;
epsilon=0.08; 
xlist_for_At=(-n_x/2+1:1:n_x/2)*dx; %xlist for A_t

%% A_T
%%% Scaling factor 2 for all cases including more than one sarcomere
M_K= exp(-((xlist_for_At-mu)/(sqrt(2)*epsilon)).^2)/(sqrt(2*pi*epsilon^2));

%% convolution 
xlist_for_concA=(0:n_x-1)*dx;
xlist_for_conv=xlist_for_At(1)+(0:2*n_x-2)*dx;
convolution_At=conv(M,M_K)*dx;

%% periodic boundary conditions 
 tleft = xlist_for_concA(1);
 tright = xlist_for_concA(end);
 jleft = find(xlist_for_conv>=tleft,1);
 jright = find(xlist_for_conv>=tright,1);
 
 xlist_for_convAt_periodic_bc=xlist_for_conv(jleft:jright); %( maybe add +/- 1) looks ok
 M_K = convolution_At(jleft:jright);
 M_K(end-(jleft-2):end) =  M_K(end-(jleft-2):end) + convolution_At(1:jleft-1);
 M_K(1:length(convolution_At)-jright) =  M_K(1:length(convolution_At)-jright) + convolution_At(jright+1:end);
end 