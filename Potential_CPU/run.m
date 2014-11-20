clear all;
clc;

Z = 79; Dim = 3; RMSPro = 0.0; sigma = sqrt(RMSPro^2/3); IntTyp = 0;

rmin = 1e-02; rmax = 5.0; nr = 64;
dlnr = log(rmax/rmin)/(nr-1); r = rmin*exp((0:1:(nr-1))*dlnr);

tic;
[V1, dV1] = getPotential(1, Z, sigma, IntTyp, Dim, r);
[V2, dV2] = getPotential(2, Z, sigma, IntTyp, Dim, r);
[V3, dV3] = getPotential(3, Z, sigma, IntTyp, Dim, r);
[V4, dV4] = getPotential(4, Z, sigma, IntTyp, Dim, r);
[V5, dV5] = getPotential(5, Z, sigma, IntTyp, Dim, r);
[V6, dV6] = getPotential(6, Z, sigma, IntTyp, Dim, r);
toc;

figure(1);
subplot(1, 2, 1);
plot(r, V3, '-b', r, V4, '-c', r, V5, '-k', r, V6, '-r');
% plot(r, V1, '-k', r, V2, '-y', r, V3, '-c', r, V4, '-b', r, V5, '-m', r, V6, '-r');
% xlim([0 rmax]);
subplot(1, 2, 2);
plot(r, dV3, '-b', r, dV4, '-c', r, dV5, '-k', r, dV6, '-r');
% % hold on
% plot(r, dV1, '-k', r, dV2, '-y', r, dV3, '-c', r, dV4, '-b', r, dV5, '-m', r, dV6, '-r');
% xlim([0 0.1]);
[V4, dV4] = getPotential(4, 14, sigma, IntTyp, 2, 0.1)