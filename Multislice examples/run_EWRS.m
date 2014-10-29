clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;            % Gpu card
TEM.SimType = 10;       % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
TEM.MulOrder = 1;       % 1: First order MS, 2: Second Order MS
TEM.nConfFP = 0;        % Number of frozen phonon configurations
TEM.DimFP = 111;        % Dimensions phonon configurations
TEM.SeedFP = 1983;       % Frozen phonon random seed
TEM.PotPar = 6;         % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;        % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;       % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;     % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;   % Zero defocus plane
TEM.ApproxModel = 1;    % 1: MS, 2: PA, 3POA, 4:WPOA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 200;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 4; na = n; nb = n; nc = 3; ncu = 16; rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3); npuc = 512;
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Au001Crystal(na, nb, nc, ncu, sigma);
% TEM.Atoms = [5 5 0.0 7 sigma 1.0; 5 10 0.0 6 sigma 1.0];
% TEM.lx = 15;  TEM.ly = 15; TEM.dz = 0.25;
TEM.nx = na*npuc; TEM.ny = nb*npuc; nxh = TEM.nx/2; nyh = TEM.ny/2;
% TEM.dz = 0.5;
tic;
[aPsi, aM2Psi] = MULTEMMat(TEM);
toc;
figure(1);
subplot(1, 2, 1);
imagesc(aM2Psi);
colormap gray;
axis image;
subplot(1, 2, 2);
imagesc(angle(aPsi));
colormap gray;
axis image;