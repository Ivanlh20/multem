clear all; clc;

global TEM;
DefaultValues;% Load default values;

TEM.gpu = 0;                % Gpu card
TEM.SimType = 10;           % 1: STEM, 2: CBED, 3: HRTEM, 4: ED, 5: PED, 6: HCI, ... 10: EW real, 11: EW Fourier
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 300;
TEM.theta = 0.0; TEM.phi = 0; % Till ilumination (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmsAu3D = 0.085; sigma = sqrt(rmsAu3D^2/3);
TEM.lx = 10;  TEM.ly = 10; TEM.dz = 0.10;
TEM.nx = 2048; TEM.ny = 2048; nxh = TEM.nx/2; nyh = TEM.ny/2;
TEM.Atoms = [TEM.lx/2 TEM.ly/2 0.00 79 sigma 1.0; TEM.lx/2 TEM.ly/2 4.00 79 sigma 1.0;
    TEM.lx/2 TEM.ly/2 8.00 79 sigma 1.0;TEM.lx/2 TEM.ly/2 12.00 79 sigma 1.0;
    TEM.lx/2 TEM.ly/2 16.00 79 sigma 1.0;TEM.lx/2 TEM.ly/2 20.00 79 sigma 1.0];
cc = [0 0 0; 0 0 1; 0 1 1; 1 0 1; 1 0 0; 0 1 0; 1 1 0]; 
     %black, blue,  cyan,  magenta, red, green
adz = [1 1/2 1/4 1/8 1/16 1/32 1/64]*4.0;
idz = 1;
figure(1); %clf;
for dz = adz;
    TEM.dz = dz;
    tic;
    clear MULTEMMat;
    [aPsi, aM2Psi] = MULTEMMat(TEM);
    toc; 
    [sum(aM2Psi(:)), sum(abs(aPsi(:)).^2)]/(TEM.nx*TEM.ny)
    hold on;
    plot(angle(aPsi(nyh+1, (nxh+1):(nxh+60))), ':', 'color', cc(idz, :));
    idz = idz + 1;
end;