function [wn] = ncon_btbc(wa,dtdx,m)
%   
%  Non conservative scheme for the Burgers equation
%
wam=zeros(1,m);
wan=zeros(1,m);
cflw=zeros(1,m);
%
wam(1:m-1)=wa(2:m);
% Transmissive boundary conditions
wam(m)=wa(m-1);
%
wan(2:m)=wa(1:m-1);
% Transmissive boundary conditions
wan(1)=wa(2);
%
cflw(1:m)=dtdx.*wa(1:m);
%
for i=1:m
    if(wa(i)>=0)
    wn(i)=wa(i)-...
          cflw(i)*( wa(i)-wan(i) );
    else
    wn(i)=wa(i)-...
          cflw(i)*( wam(i)-wa(i) );
    end
end
