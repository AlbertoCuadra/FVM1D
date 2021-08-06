function [w] = w_cauchy(x,t,P,lambda)
    u = zeros(length(lambda),1);
    for i=1:length(lambda)
        u(i,1) = u0_cauchy(x-lambda(i)*t,P,i);
    end
    w = P*u;
end
    