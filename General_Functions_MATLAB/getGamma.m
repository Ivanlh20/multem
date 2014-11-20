% Input: E0(keV), Output: gamma(relativistic factor)
function [gamma] = getGamma(E0)
emass = 510.99906;          % electron rest mass in keV
gamma = (1.0 + E0/emass);