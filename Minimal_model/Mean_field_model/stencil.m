% Laplacian is calculated using the 5-point stencil
% Copyright (C) <2023>  <Francine Kolley>
% Francine Kolley 
% Physics of Life, Benjamin M. Friedrich group
% TU_dresden 
% contact: francine.kolley@tu-dresden.de
% Latest code 07-2023
function[L]=stencil(n,dx)

z=zeros(1,n);

z(1)=-30;
z(2)=16;
z(3)=-1;
z(n)=16;
z(n-1)=-1;

L=toeplitz(z)/(12*dx^2);
end
