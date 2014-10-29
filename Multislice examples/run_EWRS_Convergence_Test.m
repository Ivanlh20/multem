clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;            % Gpu card
TEM.SimType = 10;       % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
TEM.MulOrder = 2;       % 1: First order MS, 2: Second Order MS
TEM.nConfFP = 75;        % Number of frozen phonon configurations
TEM.DimFP = 110;        % Dimensions phonon configurations
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
na = 5; nb = 5; nc = 2; ncu = 16; rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3);
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Au001Crystal(na, nb, nc, ncu, sigma);
TEM.nx = 256*na; TEM.ny = 256*nb; nxh = TEM.nx/2; nyh = TEM.ny/2;
tic;
[aPsi, aM2Psi] = MULTEMMat(TEM);
M2aPsi = abs(aPsi).^2;
aM2Psi = AverageMatrix(aM2Psi, na, nb);
M2aPsi = AverageMatrix(M2aPsi, na, nb);
[ny, nx] = size(aM2Psi); nxh = nx/2; nyh = ny/2;
toc;
figure(1);
subplot(2, 2, 1);
imagesc(aM2Psi);
colormap gray;
axis image;
subplot(2, 2, 2);
imagesc(M2aPsi);
colormap gray;
axis image;
subplot(2, 2, 3);
hold on;
plot(aM2Psi(nyh+1,:), '-g');
xlim([0 nx]);
subplot(2, 2, 4);
hold on;
plot(M2aPsi(nyh+1,:), '-g');
xlim([0 nx]);