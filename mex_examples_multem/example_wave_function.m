clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_ncores = 1; 
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;
system_conf.gpu_nstream = 1;

% eTEMST_EWFS=51, eTEMST_EWRS=52
input_multislice.simulation_type = 52;             
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.pn_dim = 110;                      % phonon dimensions
input_multislice.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 5;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.thick_type = 1;                % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multislice.thick = 0;                     % Array of thickes

input_multislice.bwl = 0;

input_multislice.E_0 = 100;
input_multislice.theta = 0.01;
input_multislice.phi = 0.0;

na = 8; nb = 8; nc = 3; ncu = 2; rms3d = 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.thick_type = 1;             % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multislice.thick = 0:2*c:1000;         % Array of thickes

input_multislice.nx = 2048; 
input_multislice.ny = 2048;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 1;       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = read_psi_0_multem(input_multislice.nx, input_multislice.ny);    % user define incident wave
input_multislice.iw_x = 0.5*input_multislice.spec_lx;          % x position 
input_multislice.iw_y = 0.5*input_multislice.spec_ly;          % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                  % Vortex momentum
input_multislice.cond_lens_c_10 = 1110;               % Defocus (Å)
input_multislice.cond_lens_c_30 = 3.3;              % Third order spherical aberration (mm)
input_multislice.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.cond_lens_c_12 = 0.0;             % Twofold astigmatism (Å)
input_multislice.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.cond_lens_c_23 = 0.0;             % Threefold astigmatism (Å)
input_multislice.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 7.50;  % Outer aperture (mrad)
input_multislice.cond_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.cond_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.cond_lens_beta = 0.2;             % Divergence semi-angle (mrad)
input_multislice.cond_lens_nbeta = 10;             % Number of integration steps for the divergence semi-angle
input_multislice.cond_lens_zero_defocus_type = 1;  % eZDT_First = 1, eZDT_User_Define = 2
input_multislice.cond_lens_zero_defocus_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.obj_lens_m = 0;                  % Vortex momentum
input_multislice.obj_lens_c_10 = 15.836;             % Defocus (Å)
input_multislice.obj_lens_c_30 = 1e-03;            % Third order spherical aberration (mm)
input_multislice.obj_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multislice.obj_lens_c_12 = 0.0;             % Twofold astigmatism (Å)
input_multislice.obj_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multislice.obj_lens_c_23 = 0.0;             % Threefold astigmatism (Å)
input_multislice.obj_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multislice.obj_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.obj_lens_outer_aper_ang = 24.0;  % Outer aperture (mrad)
input_multislice.obj_lens_sf = 32;                % Defocus Spread (Å)
input_multislice.obj_lens_nsf = 10;               % Number of integration steps for the defocus Spread
input_multislice.obj_lens_zero_defocus_type = 3;  % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.obj_lens_zero_defocus_plane = 0;

clear il_wave_function;
tic;
ouput_multislice = il_wave_function(system_conf, input_multislice);
toc;

figure(1); 
for ithk=1:length(ouput_multislice.thick)
    psi_coh = flipud(ouput_multislice.data(ithk).psi_coh);  

    subplot(1, 2, 1);    
    imagesc(abs(psi_coh).^2);
    colormap gray;
    axis image;
    title(strcat('Intensity, thick = ', num2str(ithk)));
    subplot(1, 2, 2);    
    imagesc(angle(psi_coh));
    colormap gray;
    axis image; 
    title(strcat('Phase, thick = ', num2str(ithk))); 
    pause(0.25);
end