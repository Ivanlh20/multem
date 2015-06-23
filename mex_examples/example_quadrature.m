clear all; clc;
% QuadType
% 0: int_-1^1 f(x) dx
% 1: int_0^infty f(x) dx
% 2: int_0^infty f(x)exp(-x) dx
% 3: int_-infty^infty f(x) dx
% 4: int_0^infty f(x)sin(wx) dx
% 5: int_0^infty f(x)Cos(wx) dx
% 6: int_-1^1 f(x) dx
% 7: int_0^infty f(x) Exp[-x^2] dx
% 8: int_-infty^infty f(x) Exp[-x^2] dx
% 9: int_0^infty f(x) Exp[-x] dx
% 10: int_0^infty f(x) Exp[-x]/Sqrt[x] dx

QuadType = 1; nQuad = 128;

% Load quadrature
[xi, wi] = get_quadrature(QuadType, nQuad);

Z = 6; PotPar = 6; sigma = 0.0; IntTyp = 0; Dim = 3;
[rhoi, ~] = get_Pr(6, Z, xi);
[Vi, ~] = get_Vr(PotPar, Z, xi);
% get volumen
Volume1 = sum(4*pi*(xi.^2).*Vi.*wi);

%first component of f_e(g)
cPotf=47.877645145863056; g0 = 0.0;
Volume2 = get_feg(6, Z, g0);

num2str([Volume1 cPotf*Volume2], 8)
