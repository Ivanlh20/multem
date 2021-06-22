clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_ncores = 1; % Number of Cores CPU (It will be used in the future)
system_conf.cpu_n_thread = 4; % Number of CPU threads 
system_conf.gpu_device = 0; % GPU device (i.e. 0, 1, 2, ... )

input_multem.atomic_vib_model = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.interaction_model = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_planes = 1, esst_dz_proj = 2, esst_planes_sub = 3, esst_dz_sub = 4, esst_auto = 5
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_dim = [true, true, false]; % phonon dimensions
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)
input_multem.atomic_vib_sgl_conf = 1; % 1: true, 0:false
input_multem.atomic_vib_nconf = 1; % true: phonon configuration, false: number of frozen phonon configurations

input_multem.bwl = 0;

input_multem.E_0 = 300;
input_multem.theta = 0.0;
input_multem.phi = 0.0;

na = 4;nb = 4;nc =2;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

% input_multem.spec_atoms = [79 4.0 4.0 0 rmsd_3d 1.0];
% input_multem.spec_bs_x = 8.0;
% input_multem.spec_bs_y = 8.0;
% input_multem.spec_slic(1).dz = 0.5;

input_multem.nx = 2048;
input_multem.ny = 2048;

clear ilc_spec_slicing;
[atoms, Slice] = ilc_spec_slicing(input_multem);

[natoms, ~] = size(atoms);[nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multem.islice = islice;
    
    system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 2; % eP_Float = 1, eP_double = 2
    tic;
%  clear ilc_transmission_function;
    output_multislice_1 = ilc_transmission_function(system_conf, input_multem);
    toc;
    
    system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 2; % eP_Float = 1, eP_double = 2
    tic;
%  clear ilc_transmission_function;
    output_multislice_2 = ilc_transmission_function(system_conf, input_multem);
    toc;
    sum(abs(output_multislice_1.trans(:)-output_multislice_2.trans(:))/(input_multem.nx*input_multem.ny))
    
    figure(1);
    subplot(1, 3, 1);
    imagesc(real(output_multislice_1.trans));
    colormap gray;
    axis image;
    subplot(1, 3, 2);
    imagesc(imag(output_multislice_1.trans));
    colormap gray;
    axis image;
    subplot(1, 3, 3);
    imagesc(abs(output_multislice_1.trans));
    colormap gray;
    axis image;
    num2str([islice, min(abs(output_multislice_1.trans(:))), max(abs(output_multislice_1.trans(:))), sum(abs(output_multislice_1.trans(:)))/(input_multem.nx*input_multem.ny)], 10)
    pause(0.10);
end