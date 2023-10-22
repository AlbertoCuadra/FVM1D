function [wn]=god_btbc(wa,dtdx,m)
%   
%  Godunov scheme for the Burgers equation
%
fw=zeros(1,m);
fwn=zeros(1,m);
wn=zeros(1,m);
fw(1:m)=0.5*wa(1:m).^2;
% This vaule is the correct one for flux at the left state
% and the constant if wr=wr for consistency 
  for i=1:m-1
      wl=wa(i);
      wr=wa(i+1);
      s=(wl+wr)/2;
      if wl > wr 
        % Entropy shock
        if s < 0
         fw(i) = 0.5*wr^2;
        end
      elseif wl < wr 
        % Expansion wave
        if wr < 0
         fw(i) = 0.5*wr^2 ;
        elseif wl > 0.
         fw(i) = 0.5*wl^2 ;   
        else
         fw(i) = 0;
        end
      end
  end
%
fwn(2:m)=fw(1:m-1);
% Transmissive boundary conditions
fwn(1)=fw(1);
%
wn(1:m)=wa(1:m)-dtdx*(fw(1:m)-fwn(1:m));
%
end

