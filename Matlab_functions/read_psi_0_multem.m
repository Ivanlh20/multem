function [psi] = read_psi_0_multem(nx, ny)
psi_0 = importdata('m2psi_tot_0.mat');
[ny0, nx0] = size(psi_0);
[Rx, Ry] = meshgrid(1+(0:(nx-1))*(nx0-1)/(nx-1), 1+(0:(ny-1))*(ny0-1)/(ny-1));
A = interp2(psi_0, Rx, Ry);
P = exp(1i*0.5*pi*A);
psi = A.*P;