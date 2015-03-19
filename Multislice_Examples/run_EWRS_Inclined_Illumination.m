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
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
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
na = 5; nb = 5; nc = 2; ncu = 4; rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3);
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Au001Crystal(na, nb, nc, ncu, sigma);
TEM.nx = 256*na; TEM.ny = 256*nb; nxh = TEM.nx/2; nyh = TEM.ny/2;
tic;
[aPsi, aM2Psi] = MULTEM_GPU(TEM);
figure(2);
imagesc(aM2Psi);
colormap hot;
axis image;