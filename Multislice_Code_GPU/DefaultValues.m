%Ivan lobato - 26/02/2013-11.04am
global TEM
TEM.gpu = 0;                % Gpu card
TEM.SimType = 51;           % 11: STEM, 12: ISTEM, 21: CBED, 22: CBEI, 31: ED, 32: HRTEM, 41: PED, 42: HCI, ... 51: EW Fourier, 52: EW real
TEM.nConfFP = 0;            % Number of frozen phonon configurations
TEM.DimFP = 111;            % Dimensions phonon configurations
TEM.SeedFP = 1983;          % Frozen phonon random seed
TEM.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
TEM.MEffect = 1;            % 1: Exit wave Partial coherente mode, 2: Transmission cross coefficient
TEM.STEffect = 1;           % 1: Spatial and temporal, 2: Temporal, 3: Spatial
TEM.ZeroDefTyp = 3;         % 1: First atom, 2: Middle point, 3: Last atom, 4: Fix Plane
TEM.ZeroDefPlane = 0;       % Zero defocus plane
TEM.ApproxModel = 1;        % 1: MS, 2: PA, 3:POA, 4:WPOA
TEM.BandwidthLimit = 1;     % 1: true, 2: false
TEM.ThicknessTyp = 1;       % 1: Whole specimen, 2: Throught thickness, 3: Through planes
TEM.Thickness = 0;          % Array of thicknesses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.E0 = 300;
TEM.theta = 0.0; TEM.phi = 0.0; % Till ilumination (degrees)
TEM.nx = 256; TEM.ny = 256;
TEM.lx = 10; TEM.ly = 10; 
TEM.dz = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%% Microscope effects %%%%%%%%%%%%%%%%%%%%%%%%
TEM.MC.m = 0;       %mm
TEM.MC.f = 0.0;     %Angs
TEM.MC.Cs3 = 0.00;	%mm
TEM.MC.Cs5 = 0.00;	%mm
TEM.MC.mfa2 = 0.0; TEM.MC.afa2 = 0.0; %(Angs, degrees)
TEM.MC.mfa3 = 0.0; TEM.MC.afa3 = 0.0; %(Angs, degrees)
TEM.MC.aobjl = 0.0; TEM.MC.aobju = 10000.0; %(mrad, mrad)
TEM.MC.sf = 32; TEM.MC.nsf = 10; % (Angs, number of steps)
TEM.MC.beta = 0.2; TEM.MC.nbeta = 10; %(mrad, half number of steps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.STEM.line = 0;
TEM.STEM.FastCal = 0;           % 0: normal mode(low memory consumption), 1: fast calculation(high memory consumption)
TEM.STEM.ns = 10;
TEM.STEM.x1u = 0.0; TEM.STEM.y1u = 0.0;
TEM.STEM.x2u = 1.0; TEM.STEM.y2u = 1.0;
TEM.STEM.nDet = 1;
TEM.STEM.DetCir(1).InnerAng = 60; TEM.STEM.DetCir(1).OuterAng = 180; % Inner angle(mrad) and Outer angle(mrad)
TEM.STEM.DetCir(2).InnerAng = 80; TEM.STEM.DetCir(2).OuterAng = 120; % Inner angle(mrad) and Outer angle(mrad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CBED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CBED.x0 = 0.0;      % x position 
TEM.CBED.y0 = 0.0;      % y position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CBED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.CBEI.x0 = 0.0;      % x position 
TEM.CBEI.y0 = 0.0;      % y position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HRTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.PED.nrot = 360;     % number of orientations
TEM.PED.theta = 3.0;    % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HCI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TEM.HCI.nrot = 360;         % number of orientations
TEM.HCI.theta = 3.0;        % Precession angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%% EW Real Space %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% EW Fourier Space %%%%%%%%%%%%%%%%%%%%%%
