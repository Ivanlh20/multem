clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 2;            % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
TEM.nConfFP = 0;            % Number of frozen phonon configurations
TEM.DimFP = 111;            % Dimensions phonon configurations
TEM.SeedFP = 1983;          % Frozen phonon random seed
TEM.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;            % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;           % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;       % Zero defocus plane
TEM.ApproxModel = 1;        % 1: MS, 2: PA, 3: POA, 4:WPOA
TEM.BandwidthLimit = 1;     % 1: true, 2: false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 100;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.MC.m = 0;       %mm
TEM.MC.f = 1110;     %Angs
TEM.MC.Cs3 = 3.3;	%mm
TEM.MC.Cs5 = 0.00;	%mm
TEM.MC.mfa2 = 0.0; TEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.MC.mfa3 = 0.0; TEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.MC.aobjl = 0.0; TEM.MC.aobju = 7.5; %(mrad, mrad)
TEM.MC.sf = 32; TEM.MC.nsf = 10; % (Angs, number of steps)
TEM.MC.beta = 0.2; TEM.MC.nbeta = 10; %(mrad, half number of steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
na = 8; nb = 8; nc = 55; ncu = 8; sigma = 0.076;  nxyuc = 512;
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Si001Crystal(na, nb, nc, ncu, sigma);
TEM.nx = na*nxyuc; TEM.ny = nb*nxyuc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CBED.x0 = 0.0;      % x position 
TEM.CBED.y0 = 0.0;      % y position
TEM.CBED.Space = 2;     % 1: Real space 2: Reciprocal space
tic;
[aPsi, aM2Psi] = MULTEMMat(TEM);
toc;

nxh = TEM.nx/2; d = 200; f = 1e+03;
rr = (nxh+1-d):1:(nxh+1+d);
aM2Psi = aM2Psi/max(aM2Psi(:));
M2Psis = log(1+f*aM2Psi(rr, rr));
Imaxs = max(M2Psis(:));

figure(1);
imshow(M2Psis, [0 Imaxs]);
axis square;
colormap gray;
axis image;