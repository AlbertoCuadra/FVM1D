function w0 = w0_cauchy(x)
    h0 = 1+2/5.*exp(-5.*x.^2);
    u0 = zeros(1,length(x));
%     h0 = 0.5.*exp(-(x-20).^2./8);
%     u0 = 100.*h0;
    w0 = [h0;u0];
%     w0 = [cos(x);sin(x)];
end