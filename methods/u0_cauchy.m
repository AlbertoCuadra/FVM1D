function u0 = u0_cauchy(x,P,i)
    u0 = P\w0_cauchy(x);
    u0 = u0(i,:); % In order to return the correct value
end