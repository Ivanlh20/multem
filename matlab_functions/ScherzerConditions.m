function [f, alpha] = ScherzerConditions(E_0, Cs3)
% E_0 = acceleration voltage [Kv]
Cs3 = Cs3*1e+07;
emass = 510.99906;
hc = 12.3984244;
lambda = hc/sqrt(E_0*(2*emass + E_0));
n = 3/4;
alpha = (4*(2*n-0.5)*lambda/Cs3)^(1/4); alpha = alpha*1e+03;
f = sqrt((2*n-0.5)*Cs3*lambda);