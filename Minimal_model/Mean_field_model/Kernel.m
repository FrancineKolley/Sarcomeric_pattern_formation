% Non-local interaction from titin/Sallimus - double gaussian kernel
% normalized area to 1

% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% Latest code 07-2023

function[M_T]=Kernel(l_t,M,n_x,dx,sigma)
% M = concentration of Myosin
% n_x= number of bins/grid points
% dx=bin size
% sigma=standard deviation

%create a x-vector
xlist_for_Mt=(-n_x/2+1:1:n_x/2)*dx; %xlist for M_t

%% non-local interaction M_T
M_T= (exp(-((xlist_for_Mt-l_t)/(sqrt(2)*sigma)).^2)+ exp(-((xlist_for_Mt+l_t)/(sqrt(2)*sigma)).^2))/(2*sqrt(2*pi*sigma^2));

%% convolution 
xlist_for_concA=(0:n_x-1)*dx;
xlist_for_conv=xlist_for_Mt(1)+(0:2*n_x-2)*dx;
convolution_Mt=conv(M,M_T)*dx;

%% periodic boundary conditions --> mirroring  
 tleft = xlist_for_concA(1);
 tright = xlist_for_concA(end);
 jleft = find(xlist_for_conv>=tleft,1);
 jright = find(xlist_for_conv>=tright,1);
 
 xlist_for_convAt_periodic_bc=xlist_for_conv(jleft:jright); %this is just for making sanity checks
 M_T = convolution_Mt(jleft:jright);
 M_T(end-(jleft-2):end) =  M_T(end-(jleft-2):end) + convolution_Mt(1:jleft-1);
 M_T(1:length(convolution_Mt)-jright) =  M_T(1:length(convolution_Mt)-jright) + convolution_Mt(jright+1:end);
end 