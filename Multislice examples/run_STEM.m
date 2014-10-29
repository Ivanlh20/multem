clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;            % Gpu card
TEM.SimType = 1;        % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
TEM.MulOrder = 2;       % 1: First order MS, 2: Second Order MS
TEM.nConfFP = 0;        % Number of frozen phonon configurations
TEM.DimFP = 111;        % Dimensions phonon configurations
TEM.SeedFP = 1983;      % Frozen phonon random seed
TEM.PotPar = 6;         % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;        % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;       % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;     % 1: First atom, 2: middle point, 3: last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;   % Zero defocus plane
TEM.ApproxModel = 1;    % 1: Mulstilice, 2: Projection approximation, 3: Phase object approximation, 4: Weak phase object approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 200;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
na = 15; nb = 11; nc = 4; ncu = 2; DWAu3D = 0.6373; rmsAu3D = sqrt(DWAu3D/(8*pi^2)); sigma = sqrt(rmsAu3D^2/3);
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Au110Crystal(na, nb, nc, ncu, sigma);
TEM.nx = 2560; TEM.ny = 2560; nxh = TEM.nx/2; nyh = TEM.ny/2;
figure(1);
plot3(TEM.Atoms(:, 1), TEM.Atoms(:, 2), TEM.Atoms(:, 3), '*r');
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.MC.m = 0;       %mm
TEM.MC.f = 15.836;  %Angs
TEM.MC.Cs3 = 1e-03;	%mm
TEM.MC.Cs5 = 0.00;	%mm
TEM.MC.mfa2 = 0.0; TEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.MC.mfa3 = 0.0; TEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.MC.aobjl = 0.0; TEM.MC.aobju = 24.0; %(mrad, mrad)
TEM.MC.sf = 32; TEM.MC.nsf = 10; % (Angs, number of steps)ne
TEM.MC.beta = 0.2; TEM.MC.nbeta = 10; %(mrad, half number of steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.STEM.line = 1;
TEM.STEM.FastCal = 0; 
TEM.STEM.ns = 35;
TEM.STEM.x1u = 0.0; TEM.STEM.y1u = 0.0;
TEM.STEM.x2u = a; TEM.STEM.y2u = b;
TEM.STEM.nDet = 1;
TEM.STEM.DetCir(1).InnerAng = 50; TEM.STEM.DetCir(1).OuterAng = 170; % Inner angle(mrad) and Outer angle(mrad)    
tic;
clear MULTEMMat;
STEM = MULTEMMat(TEM);
toc;
Tot = STEM.DetInt(1).Tot;
Coh = STEM.DetInt(1).Coh;
% save(strcat('Tot_', num2str(nc),'.mat'), 'Tot');
% save(strcat('Coh_', num2str(nc),'.mat'), 'Coh');    

figure(2);
subplot(2, 1, 1);
hold on;
plot(Tot, '-r');
subplot(2, 1, 2);
hold on;
plot(Coh, '-r');

% figure(2);
% subplot(2, 1, 1);
% imagesc(Tot);
% colormap gray;
% axis image;
% subplot(2, 1, 2);
% imagesc(Coh);
% colormap gray;
% axis image;