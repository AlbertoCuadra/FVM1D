function wn = god(x,A,dtdx,wa)
% [ind1,ind2] = ind2sub(size(A),find(A>0));
nx = length(x);
[P,D] = eig(A); 
A_plus = P*(0.5*(D+abs(D)))/P; % P*(0.5*(D+abs(D)))*inv(P)
A_minus= P*(0.5*(D-abs(D)))/P; % P*(0.5*(D-abs(D)))*inv(P)
% wn=zeros(length(A(:,1)),length(x));
wam(:,1:nx-1)=wa(:,2:nx);
% Transmissive boundary conditions 
wam(:,nx)=wa(:,nx-1);
wan(:,2:nx)=wa(:,1:nx-1);
% Transmissive boundary conditions
wan(:,1)=wa(:,2);
% Godunov scheme for Linearized Shallow Water 1D
wn(:,1:nx)=wa(:,1:nx)-dtdx*(A_plus*(wa(:,1:nx)-wan(:,1:nx)) + A_minus*(wam(:,1:nx)-wa(:,1:nx)));