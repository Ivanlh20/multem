clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

input_multislice.E_0 = 300;                          % Acceleration Voltage (keV)

input_multislice.lx = 10;
input_multislice.ly = 10;

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

input_multislice.lens_m = 0;           %mm
input_multislice.lens_f = 14.03;      %Angs
input_multislice.lens_Cs3 = 0.001;       %mm
input_multislice.lens_Cs5 = 0.00;       %mm
input_multislice.lens_mfa2 = 0.0; 
input_multislice.lens_afa2 = 0.0;       %(Angs, degrees)
input_multislice.lens_mfa3 = 0.0; 
input_multislice.lens_afa3 = 0.0;       %(Angs, degrees)
input_multislice.lens_aobjl = 0.0; 
input_multislice.lens_aobju = 20.0;     %(mrad, mrad)
input_multislice.lens_sf = 32; 
input_multislice.lens_nsf = 10;         % (Angs, number of steps)
input_multislice.lens_beta = 0.2; 
input_multislice.lens_nbeta = 10;       %(mrad, half number of steps)

for ap = [5 10 15 20 30 52 80 120  150]
        input_multislice.lens_aobju = ap;
        input_multislice.conv_beam_wave_x = input_multislice.lx/2;
        input_multislice.conv_beam_wave_y = input_multislice.ly/2;

        clear MULTEM;
        tic;
        output_probe = get_probe(input_multislice); 
        toc;
        probe = flipud(output_probe.probe);
        figure(2);
        subplot(1, 2, 1);
        imagesc(abs(probe).^2);
        title('intensity');
        axis image;
        colormap gray;

        subplot(1, 2, 2);
        imagesc(angle(probe));
        title('phase');
        axis image;
        colormap gray;
        pause(1.0);
end;