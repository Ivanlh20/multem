% Input: E0(keV), Output: lambda (electron wave)
function [lambda] = getLambda(E0)
emass = 510.99906;		% electron rest mass in keV
hc = 12.3984244;		% Planck's const x speed of light	
lambda = hc/sqrt(E0*(2.0*emass + E0));