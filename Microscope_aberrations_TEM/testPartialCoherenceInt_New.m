%Ivan lobato - 06/02/2013-6.30pm
clear all; clc;
global TEMIm;

TEMIm.gpu = 0;	% Gpu card
TEMIm.option = 1;   % 0: Exit wave coherente mode, 1: Transmission cross coefficient
TEMIm.Psi = importdata('D:\PHD - Simulation\Codes\Matlab\Thust\SrTiO3\sigma=0.231\Psir_pp=6_ncu=16_nc=5_k=1.mat');
TEMIm.E0 = 300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
na = 6; nb = 4; nc = 1; pp = 6; ncu = 16; rms3D = 0; sigma = sqrt(rms3D^2/3);
[Atomso, lx, ly, lz, dz, c] = SrTiO3110Crystal(na, nb, nc, ncu, sigma);
TEMIm.nx = na*256; TEMIm.ny = na*256;
TEMIm.lx = lx; TEMIm.ly = ly;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEMIm.f = 1000*10;%Angs
TEMIm.Cs3 = 0.025; %mm
TEMIm.Cs5 = 0.00; %mm
TEMIm.mfa2 = 0.0;   TEMIm.afa2 = 0.0; %(Angs, degrees)
TEMIm.mfa3 = 0.0;   TEMIm.afa3 = 0.0; %(Angs, degrees)
TEMIm.aobjl = 0.0;	TEMIm.aobju = 1000.0; % (mrad, mrad)
TEMIm.sf = 35;      TEMIm.nf = 10; % (Angs, Angs, half number of steps)
TEMIm.beta = 0.05;	TEMIm.ngs = 8; % (mrad, half number of steps)

tic;
clear PCBFTEM;
M2Psi = PCBFTEM(TEMIm);
toc;

Psi = TEMIm.Psi;
[ny, nx] = size(Psi); nxh = nx/2; dgx = 1.0/lx; nyh = ny/2; dgy = 1.0/ly;
[gx, gy] = meshgrid((-nxh:1:(nxh-1))*dgx, (-nyh:1:(nyh-1))*dgy); g2 = gx.^2+gy.^2; g2 = ifftshift(g2);
E0 = TEMIm.E0; cs = TEMIm.Cs3; f = TEMIm.f; sf = TEMIm.sf; beta = TEMIm.beta;
T = CTF(g2, 1000, E0, cs, f, sf, beta);
fPsi = fft2(ifftshift(Psi));
fPsi = ifft2(fPsi.*T); M2Psin = abs(fftshift(fPsi)).^2;  

tic;
TEMIm.option = 0; % 0: Exit wave coherente mode, 1: Transmission cross coefficient
clear PCBFTEM;
M2Psint = PCBFTEM(TEMIm);
toc;

[min(M2Psi(:)), min(M2Psin(:)), min(M2Psint(:))]
[max(M2Psi(:)), max(M2Psin(:)), max(M2Psint(:))]

sum(abs(M2Psin(:)-M2Psint(:)))
figure(2);
subplot(1, 3, 1);
Psi = AverageMatrix(Psi, na, nb);
imagesc(abs(Psi).^2);
axis square;
colormap gray;
subplot(1, 3, 2);
M2Psi = AverageMatrix(M2Psi, na, nb);
imagesc(M2Psi);
axis square;
colormap gray;
subplot(1, 3, 3);
M2Psin = AverageMatrix(M2Psin, na, nb);
imagesc(M2Psin);
axis square;
colormap gray;

M2Psin = M2Psin/mean(M2Psin(:));
M2Psi = M2Psi/mean(M2Psi(:));
[std(M2Psi(:)), std(M2Psin(:))];