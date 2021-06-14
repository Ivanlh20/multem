clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

input_multem = ilm_dflt_input_multem(); % Load default values;

system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_ncores = 1; % Number of Cores CPU (It will be used in the future)
system_conf.cpu_n_thread = 4; % Number of CPU threads 
system_conf.gpu_device = 0; % GPU device (i.e. 0, 1, 2, ... )

input_multem.E_0 = 300; % Acceleration Voltage (keV)
input_multem.theta = 0.0;
input_multem.phi = 0.0;

input_multem.nx = 1024;
input_multem.ny = 1024;

input_multem.spec_bs_x = 50;
input_multem.spec_bs_y = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 2; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = read_psi_0_multem(input_multem.nx, input_multem.ny); % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0; % Vortex momentum
input_multem.cond_lens_c_10 = -14.0312; % Defocus (Å)
input_multem.cond_lens_c_30 = 1e-03; % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.cond_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.cond_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.cond_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.cond_lens_inner_aper_ang = 0; % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0; % Outer aperture (mrad)
input_multem.cond_lens_tp_inc_sigma = 32; % standard deviation (Å)
input_multem.cond_lens_tp_inc_npts = 10; % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_spt_inc_sigma = 0.2; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.cond_lens_spt_inc_rad_npts = 8; % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_zero_defocus_type = 1; % eZDT_First = 1, eZDT_User_Define = 2
input_multem.cond_lens_zero_defocus_plane = 0;
input_multem.cond_lens_c_10 = ilc_scherzer_defocus(input_multem.E_0, input_multem.cond_lens_c_30);

input_multem.cond_lens_c_10 = input_multem.cond_lens_c_10;

tic;
output_incident_wave = ilc_incident_wave(system_conf, input_multem);
toc;

psi_0 = flipud(output_incident_wave.psi_0);
figure(1); clf;
subplot(1, 2, 1);
imagesc(output_incident_wave.x, output_incident_wave.y, abs(psi_0).^2);
title('intensity');
axis image;
colormap gray;

subplot(1, 2, 2);
imagesc(output_incident_wave.x, output_incident_wave.y, angle(psi_0));
title('phase');
axis image;
colormap gray;
pause(0.2);