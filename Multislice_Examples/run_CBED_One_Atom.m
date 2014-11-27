clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 21;           % 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
TEM.nConfFP = 0;            % Number of frozen phonon configurations
TEM.DimFP = 111;            % Dimensions phonon configurations
TEM.SeedFP = 1983;          % Frozen phonon random seed
TEM.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;            % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;           % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;       % Zero defocus plane
TEM.ApproxModel = 1;        % 1: MS, 2: PA, 3:POA, 4:WPOA
TEM.BandwidthLimit = 1;     % 1: true, 2: false
TEM.ThicknessTyp = 1;       % 1: Whole specimen, 2: Throught thickness, 3: Through planes
TEM.Thickness = 0;          % Array of thicknesses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 300;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.PBC_xy = 0;
TEM.nx = 1024; TEM.ny = 1024; nxh = TEM.nx/2; nyh = TEM.ny/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.MC.m = 0;       %mm
TEM.MC.f = 88.7414;  %Angs
TEM.MC.Cs3 = 0.04;	%mm
TEM.MC.Cs5 = 0.00;	%mm
TEM.MC.mfa2 = 0.0; TEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.MC.mfa3 = 0.0; TEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.MC.aobjl = 0.0; TEM.MC.aobju = 21.0659; %(mrad, mrad)
TEM.MC.sf = 32; TEM.MC.nsf = 10; % (Angs, number of steps)ne
TEM.MC.beta = 0.2; TEM.MC.nbeta = 10; %(mrad, half number of steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmsAu3D = 0.084; sigma = sqrt(rmsAu3D^2/3);
TEM.lx = 10; TEM.ly = 10; TEM.dz = 0.15;
TEM.Atoms = [0.5*TEM.lx, 0.5*TEM.ly, 0.0 79 sigma 1.00];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.STEM.line = 1;
TEM.STEM.ns = 100;
TEM.STEM.x1u = 0.0; TEM.STEM.y1u = 0.5*TEM.ly;
TEM.STEM.x2u = TEM.lx; TEM.STEM.y2u = 0.5*TEM.ly;
TEM.STEM.nDet = 1;
TEM.STEM.DetCir(1).InnerAng = 60; TEM.STEM.DetCir(1).OuterAng = 190; % Inner angle(mrad) and Outer angle(mrad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CBED.x0 = 0.5*TEM.lx;      % x position 
TEM.CBED.y0 = 0.5*TEM.ly;      % y position

tic;
clear MULTEMMat;
[Psi, aM2Psi] = MULTEMMat(TEM);
toc;

figure(1);
subplot(2, 1, 1);
imagesc(abs(Psi).^2);
colormap gray;
axis image;
subplot(2, 1, 2);
imagesc(abs(Psi).^2);
colormap gray;
axis image;