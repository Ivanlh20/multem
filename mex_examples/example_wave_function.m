clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

% eST_EWFS=51, eST_EWRS=52
input_multislice.simulation_type = 52;             
input_multislice.phonon_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 110;                      % phonon dimensions
input_multislice.fp_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.fp_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.zero_defocus_type = 3;             % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User = 4
input_multislice.zero_defocus_plane = 0;

input_multislice.thickness_type = 1;                % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
input_multislice.thickness = 0;                     % Array of thicknesses

input_multislice.bwl = 0;

input_multislice.E_0 = 300;
input_multislice.theta = 0.0;
input_multislice.phi = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 1;       % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define. (options 1 and 2 are only active for EWRS or EWFS)
input_multislice.iw_psi = 0;        % user define incident wave
input_multislice.iw_x = 2*4.078;    % x position 
input_multislice.iw_y = 2*4.078;    % y position

%%%%%%%%%%%%%%%%%%%%%%%%%%%% aberrations %%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.lens_m = 0;            %mm
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

na = 4; nb = 4; nc = 21; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.thickness_type = 3;             % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
input_multislice.thickness = 0:c:1000;           % Array of thicknesses

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

input_multislice.precision = 1;                  % eP_Float = 1, eP_double = 2
input_multislice.device = 2;                     % eD_CPU = 1, eD_GPU = 2
clear get_wave_function;
tic;
ouput_multislice = get_wave_function(input_multislice);
toc;

figure(1); 
for ithk=1:length(ouput_multislice.thickness)
    psi_coh = ouput_multislice.data(ithk).psi_coh;  

    subplot(1, 2, 1);    
    imagesc(abs(psi_coh).^2);
    colormap gray;
    axis image;
    title(strcat('Intensity, thickness = ', num2str(ithk)));
    subplot(1, 2, 2);    
    imagesc(angle(psi_coh));
    colormap gray;
    axis image; 
    title(strcat('Phase, thickness = ', num2str(ithk))); 
    pause(0.5);
end;