%Ivan lobato - 06/02/2013-6.30pm
clear all; clc;
global TEM;

DefaultValues;% Load default values;

TEM.gpu = 0; % Gpu card
TEM.SimType = 0; % 0: TEM, 1: STEM, 2: PED
TEM.MulOrder = 2; % Multislice order
TEM.CTEM.option = 0; % 0: Exit wave, 1: Exit wave coherente mode, 2: Transmission cross coefficient(Exit wave, intensity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.MulSli.E0 = 300;
TEM.MulSli.Rb = 0.0; TEM.MulSli.zb = 4.0; TEM.MulSli.zmax = 4.0;
TEM.MulSli.theta = 0.0; TEM.MulSli.phi = 0.0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
na = 6; nb = 4; nc = 1; pp = 6; ncu = 16; rms3D = 0; sigma = sqrt(rms3D^2/3);
[Atomso, lx, ly, lz, dz, c] = SrTiO3110Crystal(na, nb, nc, ncu, sigma);
TEM.MulSli.nx = na*256; TEM.MulSli.ny = na*256; TEM.MulSli.nz = 10;
TEM.MulSli.lx = lx; TEM.MulSli.ly = ly; TEM.MulSli.lz = lz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CTEM.option = 0; % 0: Exit wave, 1: Exit wave coherente mode, 2: Transmission cross coefficient(Exit wave, intensity)
TEM.CTEM.MC.Cs3 = -0.025; %mm
TEM.CTEM.MC.Cs5 = 0.00; %mm
TEM.CTEM.MC.mfa2 = 0.0; TEM.CTEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.CTEM.MC.mfa3 = 0.0; TEM.CTEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.CTEM.MC.aobjl = 0.0; TEM.CTEM.MC.aobju = 1000.0; % (mrad, mrad)
TEM.CTEM.MC.f = 100; TEM.CTEM.MC.sf = 33; TEM.CTEM.MC.nf = 12; % (Angs, Angs, half number of steps)
TEM.CTEM.MC.beta = 0.25; TEM.CTEM.MC.ngs = 7; % (mrad, half number of steps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi = importdata('SrTiO3\sigma =0.231\Psir_pp=6_ncu=16_nc=5_k=1.mat'); 
tic;
M2Psi = PCBFTEM(TEM, Psi);
toc;
[ny, nx] = size(Psi); nxh = nx/2; dgx = 1.0/lx; nyh = ny/2; dgy = 1.0/ly;
[gx, gy] = meshgrid((-nxh:1:(nxh-1))*dgx, (-nyh:1:(nyh-1))*dgy); g2 = gx.^2+gy.^2; g2 = ifftshift(g2);
E0 = TEM.MulSli.E0; cs = TEM.CTEM.MC.Cs3; f = TEM.CTEM.MC.f; sf = TEM.CTEM.MC.sf; beta = TEM.CTEM.MC.beta;
T = CTF(g2, 1000, E0, cs, f, sf, beta);
fPsi = fft2(ifftshift(Psi));
fPsi = ifft2(fPsi.*T); M2Psin = abs(fftshift(fPsi)).^2;  

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
[std(M2Psi(:)), std(M2Psin(:))]
% [sum(M2Psi(:)), sum(M2Psin(:))]
% sum(abs(M2Psi(:)-M2Psin(:)))
% sum(abs(M2Psic(:)-M2Psin(:)))