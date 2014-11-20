clear all;
clc;
a0 = 0.52917721077817892;
Z = 79;
rmin = 1e-05; rmax = 0.1; nr = 512;
dlnr = log(rmax/rmin)/(nr-1); r = rmin*exp((0:1:(nr-1))*dlnr);
% dr = (rmax - rmin)/(nr-1); r = rmin:dr:rmax;

[P1, dP1] = getrho(1, Z, r); 
[P2, dP2] = getrho(2, Z, r);
[P3, dP3] = getrho(3, Z, r);
[P4, dP4] = getrho(4, Z, r);
[P5, dP5] = getrho(5, Z, r);
[P6, dP6] = getrho(6, Z, r);

figure(1);
subplot(1, 2, 1);
plot(r, P3, '-b', r, P4, '-c', r, P5, '-k', r, P6, '-r');
ylim([0 max(P6)]);
subplot(1, 2, 2);
plot(r, dP3, '-b', r, dP4, '-c', r, dP5, '-k', r, dP6, '-r');
ylim([min(dP6) max(dP6)]);