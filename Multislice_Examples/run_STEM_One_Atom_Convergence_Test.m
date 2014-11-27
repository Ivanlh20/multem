clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 11;           % 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
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
TEM.MC.m = 0;       %mm
TEM.MC.f = 88.7414;  %Angs
TEM.MC.Cs3 = 0.040;	%mm
TEM.MC.Cs5 = 0.00;	%mm
TEM.MC.mfa2 = 0.0; TEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.MC.mfa3 = 0.0; TEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.MC.aobjl = 0.0; TEM.MC.aobju = 21.0659; %(mrad, mrad)
TEM.MC.sf = 32; TEM.MC.nsf = 10; % (Angs, number of steps)ne
TEM.MC.beta = 0.2; TEM.MC.nbeta = 10; %(mrad, half number of steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DWAu3D = 0.6373; rmsAu3D = sqrt(DWAu3D/(8*pi^2)); sigma = sqrt(rmsAu3D^2/3);n=2;
TEM.nx = n*1024; TEM.ny = n*1024; nxh = TEM.nx/2; nyh = TEM.ny/2;
TEM.Atoms = [5.0 5.0 0.0 79 sigma 1.0;5.0 5.0 2.5 79 sigma 1.0;5.0 5.0 5.0 79 sigma 1.0;5.0 5.0 7.5 79 sigma 1.0;5.0 5.0 10 79 sigma 1.0;5.0 5.0 12.5 79 sigma 1.0];
TEM.Atoms = [5.0 5.0 0.0 79 sigma 1.0];
TEM.lx = 10; TEM.ly = 10; TEM.dz = 0.5;

na = 20; nb = 14; nc = 10; ncu = 2; DWAu3D = 0.6373; rmsAu3D = sqrt(DWAu3D/(8*pi^2)); sigma = sqrt(rmsAu3D^2/3);
[TEM.Atoms, TEM.lx, TEM.ly, lz, a, b, c, TEM.dz] = Au110Crystal(na, nb, nc, ncu, sigma);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.STEM.line = 1;
TEM.STEM.FastCal = 1;
TEM.STEM.ns = 2;
TEM.STEM.x1u = 0.35*TEM.lx; TEM.STEM.y1u = 0.5*TEM.ly;
TEM.STEM.x2u = 0.65*TEM.lx; TEM.STEM.y2u = 0.5*TEM.ly;
% TEM.STEM.x1u = 0.0*TEM.lx; TEM.STEM.y1u = 0.0*TEM.ly;
% TEM.STEM.x2u = 1.0*TEM.lx; TEM.STEM.y2u = 1.0*TEM.ly;
TEM.STEM.nDet = 6;
TEM.STEM.DetCir(1).InnerAng = 10; TEM.STEM.DetCir(1).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(2).InnerAng = 20; TEM.STEM.DetCir(2).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(3).InnerAng = 40; TEM.STEM.DetCir(3).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(4).InnerAng = 60; TEM.STEM.DetCir(4).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(5).InnerAng = 80; TEM.STEM.DetCir(5).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(6).InnerAng = 100; TEM.STEM.DetCir(6).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(7).InnerAng = 120; TEM.STEM.DetCir(7).OuterAng = 160;	% Inner angle(mrad) and Outer angle(mrad)
% close all;
cc = [0 0 0; 0 0 1; 0 1 1; 1 0 1; 1 0 0; 0 1 0; 1 1 0]; 
     %black, blue,  cyan,  magenta, red, green
adz = [1 1/2 1/4 1/8 1/16 1/32 1/64]*2.5;
idz = 3;
% close all;
for dz = adz(2)
    TEM.dz = dz;
    tic;
    clear MULTEMMat;
    STEM = MULTEMMat(TEM);
    toc;    
    for i=1:6     
        figure(1);
        subplot(2, 3, i);
        hold on;
        plot(STEM.DetInt(i).Tot, '-','color', cc(idz, :));
%         imagesc(STEM.DetInt(i).Tot);
%         colormap gray;
    end;
    inFP = 0;
    idz = idz + 1;
%     ShowSlicing(TEM, inFP, 2);   
%     pause(0.25);
end;

% close all
