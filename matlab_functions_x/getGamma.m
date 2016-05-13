% Input: E_0(keV), Output: gamma(relativistic factor)
function [gamma] = getGamma(E_0)
emass = 510.99906;          % electron rest mass in keV
gamma = (1.0 + E_0/emass);