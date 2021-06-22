clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])


input_multem = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;

input_multem.E_0 = 200;                          % Acceleration Voltage (keV)
input_multem.theta = 0.00;
input_multem.phi = 0.0;

input_multem.spec_lx = 20;
input_multem.spec_ly = 20;

input_multem.nx = 1024; 
input_multem.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.iw_type = 2;   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.iw_psi = 0;    % user define incident wave
input_multem.iw_x = 0.0;    % x position 
input_multem.iw_y = 0.0;    % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0;                  % Vortex momentum
input_multem.cond_lens_c_10 = 0;             % Defocus (Å)
input_multem.cond_lens_c_30 = 0.002;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0;             % Twofold astigmatism (Å)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (º)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (Å)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (º)
input_multem.cond_lens_inner_aper_ang = 0.0;       % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0;      % Outer aperture (mrad)
input_multem.cond_lens_ti_sigma = 32;                % standard deviation (Å)
input_multem.cond_lens_ti_npts = 10;               % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_sigma = 0.2;             % standard deviation: For parallel ilumination(Å^-1); otherwise (Å)
input_multem.cond_lens_si_rad_npts = 8;             % # of integration points. It will be only used if illumination_model=4

input_multem.iw_x = 0.5*input_multem.spec_lx;
input_multem.iw_y = 0.5*input_multem.spec_ly;

ax = 0:input_multem.spec_lx/input_multem.nx:input_multem.spec_lx;
ay = 0:input_multem.spec_ly/input_multem.ny:input_multem.spec_ly;

df0 = ilc_scherzer_defocus(input_multem.E_0, input_multem.cond_lens_c_30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% simulation box size %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_lx = 20;
input_multem.spec_ly = 20;
input_multem.iw_x = 0.5*input_multem.spec_lx;
input_multem.iw_y = 0.5*input_multem.spec_ly;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%The incoming beam at Scherzer defocus at z=0 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_c_10 = df0;
output_incident_wave = ilc_incident_wave(system_conf, input_multem); 
psi_i = output_incident_wave.psi_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The incoming beam at Scherzer defocus after traveling a thickness(z=thk) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk = 250;
input_multem.cond_lens_c_10 = (df0+thk);
output_incident_wave = ilc_incident_wave(system_conf, input_multem); 
psi_o = output_incident_wave.psi_0;


figure(1);
subplot(1, 2, 1);
imagesc(ax, ay, abs(psi_i).^2);
title(['intensity - df = ', num2str(df0)]);
axis image;
colormap gray;
subplot(1, 2, 2);
imagesc(ax, ay, abs(psi_o).^2);
title(['intensity - df = ', num2str(df0+thk)]);
axis image;
colormap gray;
pause(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% simulation box size %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.spec_lx = 50;
input_multem.spec_ly = 50;
input_multem.iw_x = 0.5*input_multem.spec_lx;
input_multem.iw_y = 0.5*input_multem.spec_ly;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%The incoming beam at Scherzer defocus at z=0 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_c_10 = df0;
output_incident_wave = ilc_incident_wave(system_conf, input_multem); 
psi_i = output_incident_wave.psi_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The incoming beam at Scherzer defocus after traveling a thickness(z=thk) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk = 250;
input_multem.cond_lens_c_10 = (df0+thk);
output_incident_wave = ilc_incident_wave(system_conf, input_multem); 
psi_o = output_incident_wave.psi_0;


figure(2);
subplot(1, 2, 1);
imagesc(ax, ay, abs(psi_i).^2);
title(['intensity - df = ', num2str(df0)]);
axis image;
colormap gray;
subplot(1, 2, 2);
imagesc(ax, ay, abs(psi_o).^2);
title(['intensity - df = ', num2str(df0+thk)]);
axis image;
colormap gray;
pause(0.5);