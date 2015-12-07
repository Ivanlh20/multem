clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

% eST_STEM=11, eST_ISTEM=12, eST_CBED=21, eST_CBEI=22, eST_ED=31, eST_HRTEM=32, eST_PED=41, eST_HCI=42, eST_EWFS=51, eST_EWRS=52, 
% eST_EELS=61, eST_EFTEM=62, eST_ProbeFS=71, eST_ProbeRS=72, eST_PPFS=81, eST_PPRS=82,eST_TFFS=91, eST_TFRS=92
input_multislice.simulation_type = 52;
input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 110; 
input_multislice.fp_seed = 300183; 
input_multislice.fp_single_conf = 0;                % 1: true, 0:false
input_multislice.fp_nconf = 25;

input_multislice.zero_defocus_type = 3;             % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multislice.zero_defocus_plane = 0;
input_multislice.input_wave_type = 1;               % eIWT_Automatic = 1, eIWT_User_Define = 2
input_multislice.psi_0 = 0;

input_multislice.bwl = 0;

input_multislice.E_0 = 300;                         % Acceleration Voltage (keV)
input_multislice.theta = 3.0;                       % Till ilumination (degrees)
input_multislice.phi = 0.0;                         % Till ilumination (degrees)

na = 4; nb = 4; nc = 15; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Cu001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.iw_type = 1;                      % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define. (options 1 and 2 are only active for EWRS or EWFS)
input_multislice.iw_psi = 0;                       % user define incident wave
input_multislice.iw_x = input_multislice.lx/2;     % x position 
input_multislice.iw_y = input_multislice.ly/2;     % y position

%%%%%%%%%%%%%%%%%%%%%%%%%%%% lens aberrations %%%%%%%%%%%%%%%%%%%%%%%%%%%
input_multislice.lens_m = 0;        % vortex momentum
input_multislice.lens_f = 1110;     %Angs
input_multislice.lens_Cs3 = 3.3;	%mm
input_multislice.lens_Cs5 = 0.00;	%mm
input_multislice.lens_mfa2 = 0.0; 
input_multislice.lens_afa2 = 0.0; %(Angs, degrees)
input_multislice.lens_mfa3 = 0.0; 
input_multislice.lens_afa3 = 0.0; %(Angs, degrees)
input_multislice.lens_aobjl = 0.0; 
input_multislice.lens_aobju = 7.5; %(mrad, mrad)
input_multislice.lens_sf = 32; 
input_multislice.lens_nsf = 10; % (Angs, number of steps)
input_multislice.lens_beta = 0.2; 
input_multislice.lens_nbeta = 10; %(mrad, half number of steps)

input_multislice.thickness_type = 2;             % eTT_Whole_Specimen = 1, eTT_Through_Thickness = 2, eTT_Through_Slices = 3
input_multislice.thickness = 0:c:1000;           % Array of thicknesses

clear il_MULTEM;
tic;
output_multislice = il_MULTEM(input_multislice); 
toc;

figure(1);
for i=1:length(output_multislice.data)
    m2psi_coh = abs(output_multislice.data(i).psi_coh).^2;
    if(isfield(output_multislice.data, 'm2psi_tot'))
        m2psi_tot = output_multislice.data(i).m2psi_tot;  
    else
        m2psi_tot = m2psi_coh;
    end;    
    subplot(1, 2, 1);
    imagesc(m2psi_tot);
    title(strcat('Total intensity -  Thickness = ', num2str(i)));
    axis image;
    colormap gray;

    subplot(1, 2, 2);
    imagesc(m2psi_coh);
    title(strcat('Coherent intensity -  Thickness = ', num2str(i)));
    axis image;
    colormap gray;
    pause(0.25);
end;