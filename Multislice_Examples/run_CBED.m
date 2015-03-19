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
TEM.ApproxModel = 2;        % 1: MS, 2: PA, 3:POA, 4:WPOA
TEM.BWL = 1;                % 1: true, 2: false
TEM.FastCal = 1;            % 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)
TEM.ThkTyp = 1;             % 1: Whole specimen, 2: Throught thickness, 3: Through planes
TEM.Thk = 0;                % Array of thicknesses
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
na = 5; nb = 5; nc = 2; ncu = 4; sigma = 0.076;  nxyuc = 512;
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Si001Crystal(na, nb, nc, ncu, sigma);
Atoms = TEM.Atoms;
Atoms = [TEM.lx, TEM.ly, TEM.dz, 0, 0, 0; Atoms];
% save('Si001_5x5x5.txt','Atoms','-ascii');
TEM.nx = 2048; TEM.ny = 2048;
num2str([TEM.lx, TEM.ly], 8)
% na = 8; nb = 8; nc = 10; ncu = 8; sigma = 0.076;  nxyuc = 512;
% [TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Si001Crystal(na, nb, nc, ncu, sigma);
% TEM.nx = na*nxyuc; TEM.ny = nb*nxyuc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CBED.x0 = 0.0;      % x position 
TEM.CBED.y0 = 0.0;      % y position
tic;
[aPsi, aM2Psi] = MULTEM_GPU(TEM);
toc;
a = aM2Psi((1024-5):(1024+5),(1024-5):(1024+5));
% num2str(a, 10)

nxh = TEM.nx/2; d = 200; f = 1e+05;
rr = (nxh+1-d):1:(nxh+1+d);
aM2Psi = aM2Psi/max(aM2Psi(:));
M2Psis = log(1+f*aM2Psi(rr, rr));
Imaxs = max(M2Psis(:));

figure(2);
imshow(M2Psis, [0 Imaxs]);
% imshow(aM2Psi);
axis square;
colormap hot;
axis image;

% a = aM2Psi((1024-5):(1024+5),(1024-5):(1024+5));
% num2str(a, 6)