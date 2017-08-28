clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;

input_multislice.E_0 = 200;                          % Acceleration Voltage (keV)
input_multislice.theta = 0.0;
input_multislice.phi = 0.0;

input_multislice.spec_lx = 50;
input_multislice.spec_ly = 50;

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 3;   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multislice.iw_psi = read_psi_0_multem(input_multislice.nx, input_multislice.ny);    % user define incident wave
input_multislice.iw_x = 0.0;    % x position 
input_multislice.iw_y = 0.0;    % y position

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.cond_lens_m = 0;                  % Vortex momentum
input_multislice.cond_lens_c_10 = -15.836;         % Defocus (Å)
input_multislice.cond_lens_c_30 = 1e-03;           % Third order spherical aberration (mm)
input_multislice.cond_lens_c_50 = 0.00;            % Fifth order spherical aberration (mm)
input_multislice.cond_lens_c_12 = 0.0;             % Twofold astigmatism (Å)
input_multislice.cond_lens_phi_12 = 0.0;           % Azimuthal angle of the twofold astigmatism (º)
input_multislice.cond_lens_c_23 = 0.0;             % Threefold astigmatism (Å)
input_multislice.cond_lens_phi_23 = 0.0;           % Azimuthal angle of the threefold astigmatism (º)
input_multislice.cond_lens_inner_aper_ang = 0.0;   % Inner aperture (mrad) 
input_multislice.cond_lens_outer_aper_ang = 21.0;  % Outer aperture (mrad)
input_multislice.cond_lens_dsf_sigma = 32;                % standard deviation (Å)
input_multislice.cond_lens_dsf_npoints = 10;               % # of integration points. It will be only used if illumination_model=4
input_multislice.cond_lens_ssf_sigma = 0.2;             % standard deviation: For parallel ilumination(Å^-1); otherwise (Å)
input_multislice.cond_lens_ssf_npoints = 8;             % # of integration points. It will be only used if illumination_model=4

for x = (0.4:0.025:0.6)*input_multislice.spec_lx
    for y = (0.4:0.025:0.6)*input_multislice.spec_ly
        
        input_multislice.iw_x = x;
        input_multislice.iw_y = y;
%         input_multislice.iw_x = input_multislice.spec_lx/2;
%         input_multislice.iw_y = input_multislice.spec_ly/2;
        
        tic;
        output_incident_wave = il_incident_wave(system_conf, input_multislice); 
        toc;
        psi_0 = flipud(output_incident_wave.psi_0);
        figure(2);
        subplot(1, 2, 1);
        imagesc(abs(psi_0).^2);
        title('intensity');
        axis image;
        colormap gray;

        subplot(1, 2, 2);
        imagesc(angle(psi_0));
        title('phase');
        axis image;
        colormap gray;
        pause(0.1);
    end
end