clc; clear all;

clear all; clc;

Z = 79; charge = 0;
Rmin = 1e-02; Rmax = 5.0; nR = 64;
dlnR = log(Rmax/Rmin)/(nR-1); 
R = Rmin*exp((0:1:(nR-1))*dlnR);
tic;
[VR, ~] = il_Vp(1, Z, charge, R);

nR = 128;
dR = (Rmax-Rmin)/nR; 
Ri = Rmin:dR:Rmax;

VRi = il_cubic_spline(R, VR, Ri);
VRin = spline(R, VR, Ri);

figure(1);
plot(R, VR, '-+k', Ri, VRi, '-+r', Ri, VRin, '-+b');