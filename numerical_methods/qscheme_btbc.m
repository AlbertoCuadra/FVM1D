function [wn] = qscheme_btbc(wa,dtdx,m)
%   
%  Q-schme scheme for the Burgers equation
%
wam=zeros(1,m);
wan=zeros(1,m);
dm=zeros(1,m);
dn=zeros(1,m);
%
wam(1:m-1)=wa(2:m);
% Transmissive boundary conditions 
wam(m)=wa(m-1);
%
wan(2:m)=wa(1:m-1);
% Transmissive boundary conditions
wan(1)=wa(2);
%
% Central flux
%
cf(1:m)=0.25*dtdx*(wam(1:m).^2-wan(1:m).^2);
%
% Numerical viscosity Q-scheme
%
dm(1:m)=(0.25*dtdx)*abs(wa(1:m)+wam(1:m)).*(wam(1:m)-wa(1:m));
dn(1:m)=(0.25*dtdx)*abs(wa(1:m)+wan(1:m)).*(wa(1:m)-wan(1:m));
%
wn(1:m)=wa(1:m)-cf(1:m)+dm(1:m)-dn(1:m);
end

