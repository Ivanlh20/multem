% output_multislice = input_multem.ilc_multem perform TEM simulation
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])


input_multem = multem_input.parameters;         % Load default values;

input_multem.system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 4; 
input_multem.system_conf.gpu_device = 0;

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
input_multem.cond_lens_c_10 = 0;             % Defocus (�)
input_multem.cond_lens_c_30 = 0.002;            % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00;             % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0;             % Twofold astigmatism (�)
input_multem.cond_lens_phi_12 = 0.0;             % Azimuthal angle of the twofold astigmatism (�)
input_multem.cond_lens_c_23 = 0.0;             % Threefold astigmatism (�)
input_multem.cond_lens_phi_23 = 0.0;             % Azimuthal angle of the threefold astigmatism (�)
input_multem.cond_lens_inner_aper_ang = 0.0;       % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 21.0;      % Outer aperture (mrad)
input_multem.cond_lens_ti_sigma = 32;                % standard deviation (�)
input_multem.cond_lens_ti_npts = 10;               % # of integration points. It will be only used if illumination_model=4
input_multem.cond_lens_si_sigma = 0.2;             % standard deviation: For parallel ilumination(�^-1); otherwise (�)
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
output_incident_wave = input_multem.ilc_incident_wave; 
psi_i = output_incident_wave.psi_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The incoming beam at Scherzer defocus after traveling a thickness(z=thk) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk = 250;
input_multem.cond_lens_c_10 = (df0+thk);
output_incident_wave = input_multem.ilc_incident_wave; 
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
output_incident_wave = input_multem.ilc_incident_wave; 
psi_i = output_incident_wave.psi_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The incoming beam at Scherzer defocus after traveling a thickness(z=thk) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thk = 250;
input_multem.cond_lens_c_10 = (df0+thk);
output_incident_wave = input_multem.ilc_incident_wave; 
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