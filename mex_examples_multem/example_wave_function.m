clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])


input_multem = ilm_dflt_input_multem(); % Load default values;

system_config.precision = 1; % eP_Float = 1, eP_double = 2
system_config.device = 2; % eD_CPU = 1, eD_GPU = 2
system_config.cpu_n_thread = 4;
system_config.gpu_device = 0;

% eTEMST_EWFS=51, eTEMST_EWRS=52
input_multem.em_sim_typ = 52;
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6

input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)
input_multem.atomic_vib_sgl_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.atomic_vib_nconf = 5; % true: phonon configuration, false: number of frozen phonon configurations

input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0; % Array of thickes

input_multem.bwl = 0;

input_multem.E_0 = 100;
input_multem.theta = 0.01;
input_multem.phi = 0.0;

na = 8;nb = 8;nc = 3;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.thick_typ = 1; % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
input_multem.thick = 0:2*c:1000; % Array of thickes

input_multem.nx = 2048;
input_multem.ny = 2048;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.incdt_wav_typ = 1; % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto
input_multem.incdt_wav_psi = read_psi_0_multem(input_multem.nx, input_multem.ny); % user define incident wave

%%%%%%%%%%%%%%%%%%%%%%%%%%% beam position %%%%%%%%%%%%%%%%%%%%%%%%%%
input_multem.beam_pos = [input_multem.spec_bs_x/2;input_multem.spec_bs_y/2];

%%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.cond_lens_m = 0; % Vortex momentum
input_multem.cond_lens_c_10 = 1110; % Defocus (Å)
input_multem.cond_lens_c_30 = 3.3; % Third order spherical aberration (mm)
input_multem.cond_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.cond_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.cond_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.cond_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.cond_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.cond_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.cond_lens_outer_aper_ang = 7.50; % Outer aperture (mrad)
input_multem.cond_lens_tp_inc_sigma = 32; % standard deviation (Å)
input_multem.cond_lens_tp_inc_npts = 10; % # of integration points. It will be only used if illum_mod=4
input_multem.cond_lens_spt_inc_sigma = 0.2; % standard deviation: For parallel ilumination(Å^-1);otherwise (Å)
input_multem.cond_lens_spt_inc_rad_npts = 8; % # of integration points. It will be only used if illum_mod=4
input_multem.cond_lens_zero_def_typ = 1; % eZDT_First = 1, eZDT_User_Define = 2
input_multem.cond_lens_zero_def_plane = 0;

%%%%%%%%%%%%%%%%%%%%%%%% Objective lens %%%%%%%%%%%%%%%%%%%%%%%%
input_multem.obj_lens_m = 0; % Vortex momentum
input_multem.obj_lens_c_10 = 15.836; % Defocus (Å)
input_multem.obj_lens_c_30 = 1e-03; % Third order spherical aberration (mm)
input_multem.obj_lens_c_50 = 0.00; % Fifth order spherical aberration (mm)
input_multem.obj_lens_c_12 = 0.0; % Twofold astigmatism (Å)
input_multem.obj_lens_phi_12 = 0.0; % Azimuthal angle of the twofold astigmatism (º)
input_multem.obj_lens_c_23 = 0.0; % Threefold astigmatism (Å)
input_multem.obj_lens_phi_23 = 0.0; % Azimuthal angle of the threefold astigmatism (º)
input_multem.obj_lens_inner_aper_ang = 0.0; % Inner aperture (mrad) 
input_multem.obj_lens_outer_aper_ang = 24.0; % Outer aperture (mrad)
input_multem.obj_lens_tp_inc_sigma = 32; % standard deviation (Å)
input_multem.obj_lens_tp_inc_npts = 10; % # of integration points. It will be only used if illum_mod=4
input_multem.obj_lens_zero_def_typ = 3; % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4
input_multem.obj_lens_zero_def_plane = 0;

input_multem.output_area_ip_0 = [1;1]; % Starting position in pixels
input_multem.output_area_ip_e = [1;1]; % End position in pixels

clear ilc_wave_function;
tic;
output_multislice = ilc_wave_function(system_config, input_multem);
toc;

figure(1);
for ithk=1:length(output_multislice.thick)
    psi_coh = flipud(output_multislice.data(ithk).psi_coh);

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