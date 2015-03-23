clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 52;           % 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
TEM.nConfFP = 0;            % Number of frozen phonon configurations
TEM.DimFP = 111;            % Dimensions phonon configurations
TEM.SeedFP = 1983;          % Frozen phonon random seed
TEM.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;            % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;           % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;       % Zero defocus plane
TEM.ApproxModel = 1;        % 1: MS, 2: PA, 3:POA, 4:WPOA
TEM.BWL = 1;                % 1: true, 2: false
TEM.FastCal = 1;            % 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)
TEM.ThkTyp = 1;             % 1: Whole specimen, 2: Throught thickness, 3: Through planes
TEM.Thk = 0;                % Array of thicknesses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 200;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2; na = n; nb = n; nc = 3; ncu = 8; rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3); npuc = 512;
[TEM.Atoms, TEM.lx, TEM.ly, lz, TEM.dz, c] = Au001Crystal(na, nb, nc, ncu, sigma);
TEM.nx = na*npuc; TEM.ny = nb*npuc; nxh = TEM.nx/2; nyh = TEM.ny/2;
TEM.dz = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[aPsi, aM2Psi] = MULTEMMat(TEM);
zProp = -0.5*lz;
tic;
Psip = get_Propagate_GPU(aPsi, TEM.gpu, TEM.E0, TEM.nx, TEM.ny, TEM.lx, TEM.ly, zProp, TEM.BWL);
toc;
M2PsipGPU =  abs(Psip).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, g2, gmax] = g_Grid(TEM.nx, TEM.ny, TEM.lx, TEM.ly, 1);
P = Propagator(g2, gmax, TEM.E0, zProp);
fPsi = fft2(fftshift(aPsi));
fPsi = fftshift(ifft2(P.*fPsi));
M2PsipCPU = abs(fPsi).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1, 3, 1);
imagesc(aM2Psi);
colormap gray;
axis image;
TEM.ZeroDefTyp = 2;     % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
[aPsi, aM2Psi] = MULTEM_GPU(TEM);
subplot(1, 3, 2);
imagesc(M2PsipGPU);
colormap gray;
axis image;
subplot(1, 3, 3);
imagesc(M2PsipCPU);
colormap gray;
axis image;
sum(abs(M2PsipGPU(:)-M2PsipCPU(:)))

