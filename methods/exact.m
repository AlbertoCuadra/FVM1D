function [w] = exact(x,t,xb,tb,wl,wr,s)
for i=1:length(x)
    if (x(i)<=t && t<tb)
        w(i) = wl;
    elseif (t<x(i) && x(i)<=xb)
        w(i) = (wl-x(i))/(tb-t);
    elseif (x(i)>1 && t<=tb)
        w(i) = wr;
    % Entropic solution for Burgers equation
    elseif (t>=tb)
        % Entropy shock wl>wr
        if wl >= wr
            if (x(i)-xb)/(t-tb) <= s
                w(i) = wl; 
            else
                w(i) = wr;
            end
        else
        % Expansion wave wl<wr
            if (x(i)-xb)/(t-tb) <wl
                w(i) = wl;
            elseif (x(i)-xb)/(t-tb) >=wr
                w(i) = wr;          
            else
                w(i) = (x(i)-xb)/(t-tb);
            end
        end
    end
end