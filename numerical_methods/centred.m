function [wn] = centred(wa,dtdx,m)
wam(1:m-1)=wa(2:m);
wan(2:m)=wa(1:m-1);
%
% wn(2:m)=wa(2:m)-dtdx*0.5*((wam(2:m)-wan(2:m))/2).^2;
wn(2:m-1)=wa(2:m-1)-dtdx*0.5*wa(2:m-1).*((wam(2:m-1)-wan(2:m-1)));
% Transmissive boundary conditions
wn(1)=wn(2);
wn(m)=wn(m-1);
end
