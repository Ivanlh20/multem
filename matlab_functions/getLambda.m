% Input: E_0(keV), Output: lambda (electron wave)
function [lambda] = getLambda(E_0)
emass = 510.99906;		% electron rest mass in keV
hc = 12.3984244;		% Planck's const x speed of light	
lambda = hc/sqrt(E_0*(2.0*emass + E_0));