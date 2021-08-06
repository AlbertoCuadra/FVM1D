function [wn] = qscheme_llfr_btbc(wa,dtdx,m)
%   
%  Q-scheme for the Burgers equation
%  Local Lax Friedrish (LLF) regulatization
%
wam=zeros(1,m);
wan=zeros(1,m);
dm=zeros(1,m);
dn=zeros(1,m);
regm=zeros(1,m);
regn=zeros(1,m);
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
cf(1:m)=0.25*dtdx*...
          (wam(1:m).^2-wan(1:m).^2);
%
% Numerical viscosity Q-scheme with Harten regularization
%
regm(1:m)=max(abs(wa(1:m)),abs(wam(1:m)));
regn(1:m)=max(abs(wa(1:m)),abs(wan(1:m)));
%
dm(1:m)=0.5*dtdx*regm.*(wam(1:m)-wa(1:m));
%
dn(1:m)=0.5*dtdx*regn.*(wa(1:m)-wan(1:m));
%
wn(1:m)=wa(1:m)-cf(1:m)+dm(1:m)-dn(1:m);
end
