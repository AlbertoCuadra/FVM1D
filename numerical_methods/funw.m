function [w] = funw(x,t,lambda,wl,wr,alphal,alphar,P)
w = zeros(length(x),length(lambda));
for i=1:length(x)
    if (x(i)-lambda(1)*t)<=0% Case 1:
        w(i,:) = wl;
    elseif (x(i)-lambda(1)*t)*(x(i)-lambda(2)*t)<0 % Case 2:
        w(i,:) = wl + (alphar(1)-alphal(1))*P(:,1);
    elseif (x(i)-lambda(2)*t)>=0 % Case 3:
        w(i,:) = wr;
    end
end
% CHECK
% wll = 0;
% for i=1:2
% wll = wll + alphal(i)*P(:,i);
% end