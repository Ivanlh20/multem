% Input: E0(keV), Output: sigma (Interaction parameter)
function [sigma] = getSigma(E0)
emass = 510.99906;          % electron rest mass in keV
hc = 12.3984244;			% Planck's const x speed of light
x = (emass + E0)/(2.0*emass + E0);	
lambda = hc/sqrt(E0*(2.0*emass + E0));
sigma = 2.0*pi*x/(lambda*E0);