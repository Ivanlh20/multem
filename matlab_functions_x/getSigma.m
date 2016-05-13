% Input: E_0(keV), Output: sigma (Interaction parameter)
function [sigma] = getSigma(E_0)
emass = 510.99906;          % electron rest mass in keV
hc = 12.3984244;			% Planck's const x speed of light
x = (emass + E_0)/(2.0*emass + E_0);	
lambda = hc/sqrt(E_0*(2.0*emass + E_0));
sigma = 2.0*pi*x/(lambda*E_0);