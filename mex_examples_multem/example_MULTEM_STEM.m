% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Scanning transmission electron microscopy (STEM) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_ncores = 1; 
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;
system_conf.gpu_nstream = 1;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCI=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multislice.simulation_type = 11;
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.pn_dim = 110;                      % phonon dimensions (xyz)
input_multislice.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 10;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 0.0;                       % Till ilumination (º)
input_multislice.phi = 0.0;                         % Till ilumination (º)

na = 8; nb = 8; nc = 25; ncu = 2; rms3d = 0.0942; % 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 4096; 
input_multislice.ny = 4096;

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                  % Vortex momentum
input_multislice.cond_lens_c_10 = 14.0312;            % Defocus (Å)
input_multislice.cond_lens_c_30 = 1e-03;            % Third order spherical aberration (mm)
input_multislice.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.cond_lens_c_12 = 0.0;             % Twofold astigmatism (Å)
input_multislice.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.cond_lens_c_23 = 0.0;             % Threefold astigmatism (Å)
input_multislice.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)
input_multislice.cond_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.cond_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.cond_lens_beta = 0.2;             % Divergence semi-angle (mrad)
input_multislice.cond_lens_nbeta = 10;             % Number of integration steps for the divergence semi-angle
input_multislice.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.scanning_type = 2; % eST_Line = 1, eST_Area = 2
input_multislice.scanning_periodic = 1;     % 1: true, 0:false (periodic boundary conditions)
input_multislice.scanning_ns = 20;
input_multislice.scanning_x0 = 3*a; 
input_multislice.scanning_y0 = 3*b;
input_multislice.scanning_xe = 4*a;
input_multislice.scanning_ye = 4*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.detector.type = 1;  % eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
input_multislice.detector.cir(1).inner_ang = 40;  % Inner angle(mrad) 
input_multislice.detector.cir(1).outer_ang = 160; % Outer angle(mrad)
input_multislice.detector.cir(2).inner_ang = 80;  % Inner angle(mrad) 
input_multislice.detector.cir(2).outer_ang = 160; % Outer angle(mrad)

input_multislice.detector.radial(1).x  = 0;          % radial detector angle(mrad)
input_multislice.detector.radial(1).fx = 0;         % radial sensitivity value
input_multislice.detector.matrix(1).R  = 0;          % 2D detector angle(mrad)
input_multislice.detector.matrix(1).fR = 0;         % 2D sensitivity value

input_multislice.thick_type = 2;             % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multislice.thick = c:c:1000;           % Array of thickes

clear il_MULTEM;
tic;
output_multislice = il_MULTEM(system_conf, input_multislice); 
toc;
clear il_MULTEM;

figure(1);
for i=1:length(output_multislice.data)
	ndet = length(input_multislice.detector.cir);
    for j=1:ndet
        subplot(1, ndet, j);
        imagesc(output_multislice.data(i).image_tot(j).image);
        title(strcat('Thk = ', num2str(i), ', det = ', num2str(j)));
        axis image;
        colormap jet;
    end
    pause(0.25);
end
