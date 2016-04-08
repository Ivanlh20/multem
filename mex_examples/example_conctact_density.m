clear all;
clc;
a0 = 0.52917721077817892;
Z = 29; 
charge = 0;
r = 0;

[P5, dP5] = il_Pr(5, Z, charge, r);
[P6, dP6] = il_Pr(6, Z, charge, r);

t = [P5(1), P6(1)];
num2str(t, 10)
t = [P5(1), P6(1)]*a0^3;
num2str(t, 10)
P6(1)/P5(1)