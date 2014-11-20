clear all;
clc;
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

QuadType = 1; nQuad = 4;
[xi, wi] = getQuadrature(QuadType, nQuad);
% num2str([xi, wi], 20)
a0 = 0.52917721077817892;

for Z = 1:54
    [rhoi, ~] = getrho(6, Z, xi);
    ar2(Z) = sum(4*pi*(xi.^4).*rhoi.*wi/a0.^2);
end;
num2str(ar2', 6)