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
TEM.MulSli.nx = 1024; TEM.MulSli.ny = 1024; TEM.MulSli.nz = 10;
TEM.MulSli.lx = 11.7150; TEM.MulSli.ly = 16.5675; TEM.MulSli.lz = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CTEM.option = 0; % 0: Exit wave, 1: Exit wave coherente mode, 2: Transmission cross coefficient(Exit wave, intensity)
TEM.CTEM.MC.Cs3 = 0.04; %mm
TEM.CTEM.MC.Cs5 = 0.00; %mm
TEM.CTEM.MC.mfa2 = 0.0; TEM.CTEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.CTEM.MC.mfa3 = 0.0; TEM.CTEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.CTEM.MC.aobjl = 0.0; TEM.CTEM.MC.aobju = 1000.0; % (mrad, mrad)
TEM.CTEM.MC.f = 0; TEM.CTEM.MC.sf = 10; TEM.CTEM.MC.nf = 16; % (Angs, Angs, half number of steps)
TEM.CTEM.MC.beta = 0.25; TEM.CTEM.MC.ngs = 10; % (mrad, half number of steps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi = importdata('Psi-13.mat'); 
tic;
[M2Psi, M2Psic] = PCInt(TEM, Psi);
% [M2Psi] = PCInt(TEM, Psi);
toc;
% save('M2Psi-300.mat', 'M2Psi'); 
[ny, nx] = size(Psi);
nxh = nx/2; dgx = 1.0/TEM.MulSli.lx; nyh = ny/2; dgy = 1.0/TEM.MulSli.ly;
[gx, gy] = meshgrid((-nxh:1:(nxh-1))*dgx, (-nyh:1:(nyh-1))*dgy); g2 = gx.^2+gy.^2; g2 = ifftshift(g2);
E0 = TEM.MulSli.E0; cs = TEM.CTEM.MC.Cs3; f = TEM.CTEM.MC.f; sf = TEM.CTEM.MC.sf; beta = TEM.CTEM.MC.beta;
T = CTF(g2, E0, cs, f, sf, beta);
fPsi = fft2(ifftshift(Psi));
fPsi = ifft2(fPsi.*T); M2Psin = abs(fftshift(fPsi)).^2;  
figure(2);
subplot(2, 2, 1);
imagesc(abs(Psi).^2);
axis square;
colormap gray;
subplot(2, 2, 2);
imagesc(M2Psi);
axis square;
colormap gray;
subplot(2, 2, 3);
imagesc(M2Psin);
axis square;
colormap gray;
subplot(2, 2, 4);
imagesc(M2Psic);
axis square;
colormap gray;

M2Psin = M2Psin/mean(M2Psin(:));
std(M2Psin(:))
M2Psic = M2Psic/mean(M2Psic(:));
std(M2Psic(:))
% sum(abs(M2Psi(:)-M2Psin(:)))
% sum(abs(M2Psic(:)-M2Psin(:)))