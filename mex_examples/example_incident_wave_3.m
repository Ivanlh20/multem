clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

input_multislice.E_0 = 300;                          % Acceleration Voltage (keV)
input_multislice.theta = 0.0;
input_multislice.phi = 0.0;

input_multislice.lx = 8*4.079;
input_multislice.ly = 8*4.079;

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 2;   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define(options 1 and 2 are only active for EWRS or EWFS)
input_multislice.iw_psi = read_psi_0_multem(input_multislice.nx, input_multislice.ny);    % user define incident wave
input_multislice.iw_x = 0.0;    % x position 
input_multislice.iw_y = 0.0;    % y position

%%%%%%%%%%%%%%%%%%%%%%%% Microscope effects %%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.lens_m = 0;       % vortex momentum
input_multislice.lens_f = 15.836;       %Angs
input_multislice.lens_Cs3 = 1e-03;      %mm
input_multislice.lens_Cs5 = 0.00;       %mm
input_multislice.lens_mfa2 = 0.0; 
input_multislice.lens_afa2 = 0.0;       %(Angs, degrees)
input_multislice.lens_mfa3 = 0.0; 
input_multislice.lens_afa3 = 0.0;       %(Angs, degrees)
input_multislice.lens_aobjl = 0.0; 
input_multislice.lens_aobju = 24.0;     %(mrad, mrad)
input_multislice.lens_sf = 32; 
input_multislice.lens_nsf = 10;         % (Angs, number of steps)
input_multislice.lens_beta = 0.2; 
input_multislice.lens_nbeta = 10;       %(mrad, half number of steps)

input_multislice.iw_x = 0.5*input_multislice.lx;
input_multislice.iw_y = 0.5*input_multislice.ly;

df0 = 15.836;
for df = (df0):8:(df0+40*4.079)
    input_multislice.lens_f = df;      %Angs
    clear il_MULTEM;
    tic;
    output_incident_wave = il_incident_wave(input_multislice); 
    toc;
    psi_0 = flipud(output_incident_wave.psi_0);

    figure(2);
    subplot(1, 2, 1);
    imagesc(abs(psi_0).^2);
    title(strcat('intensity - df = ', num2str(df)));
    axis image;
    colormap gray;

    subplot(1, 2, 2);
    imagesc(angle(psi_0));
    title(strcat('phase - df = ', num2str(df)));
    axis image;
    colormap gray;
    pause(0.25);
end;