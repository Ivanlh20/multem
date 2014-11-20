clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 10;           % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
TEM.MulOrder = 1;           % 1: First order MS, 2: Second Order MS
TEM.nConfFP = 0;            % Number of frozen phonon configurations
TEM.DimFP = 110;            % Dimensions phonon configurations
TEM.SeedFP = 1983;          % Frozen phonon random seed
TEM.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;            % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;           % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;       % Zero defocus plane
TEM.ApproxModel = 2;        % 1: MS, 2: PA, 3:POA, 4:WPOA
TEM.BandwidthLimit = 1;     % 1: true, 2: false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 200;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2; na = n; nb = n; nc = 3; ncu = 8; rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3); npuc = 256;
[TEM.Atoms, TEM.lx, TEM.ly, lz, TEM.dz, c] = Au001Crystal(na, nb, nc, ncu, sigma);
TEM.nx = na*npuc; TEM.ny = nb*npuc; nxh = TEM.nx/2; nyh = TEM.ny/2;
TEM.dz = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Psi, aM2Psi] = MULTEMMat(TEM);
zProp = 4;
tic;
PsipGPU = getPropagate(Psi, TEM.gpu, TEM.E0, TEM.nx, TEM.ny, TEM.lx, TEM.ly, zProp, TEM.BandwidthLimit);
toc;
M2PsipGPU =  abs(PsipGPU).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, g2, gmax] = g_Grid(TEM.nx, TEM.ny, TEM.lx, TEM.ly, 1);
P = Propagator(g2, gmax, TEM.E0, zProp);
fPsi = fft2(fftshift(Psi));
fPsi = fftshift(ifft2(P.*fPsi));
M2PsipCPU = abs(fPsi).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1, 3, 1);
imagesc(abs(Psi).^2);
colormap gray;
axis image;
subplot(1, 3, 2);
imagesc(M2PsipGPU);
colormap gray;
axis image;
subplot(1, 3, 3);
imagesc(M2PsipCPU);
colormap gray;
axis image;
sum(abs(M2PsipGPU(:)-M2PsipCPU(:)))

