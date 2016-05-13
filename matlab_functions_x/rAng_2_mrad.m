% Input: E_0(keV), Output: theta in reciprocal Angstrom
function [theta] = rAng_2_mrad(rAng, E_0)
emass = 510.99906;		% electron rest mass in keV
hc = 12.3984244;		% Planck's const x speed of light	
lambda = hc/sqrt(E_0*(2.0*emass + E_0));
theta = rAng*lambda/1e-03;