clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_n_thread = 4;
system_conf.gpu_device = 0;

input_multem.atomic_vib_mod = 3; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interac_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 2; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_pot_parm_typ = 6; % eappt_doyle_0_4 = 1, eappt_peng_0_4 = 2, eappt_peng_0_12 = 3, eappt_kirkland_0_12 = 4, eappt_weickenmeier_0_12 = 5, eappt_lobato_0_12 = 6

input_multem.atomic_vib_dim = 111; % phonon dimensions
input_multem.atomic_vib_seed = 300183; % Random seed(frozen phonon)
input_multem.atomic_vib_sgl_conf = 1; % 1: true, 0:false
input_multem.atomic_vib_nconf = 1; % true: phonon configuration, false: number of frozen phonon configurations

na = 4;nb = 4;nc = 4;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.nx = 2048;
input_multem.ny = 2048;

clear ilc_spec_slicing;
[atoms, Slice] = ilc_spec_slicing(input_multem);

[natoms, ~] = size(atoms);[nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multem.islice = islice;
    
    system_conf.device = 1; % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 1; % eP_Float = 1, eP_double = 2
    tic;
    clear ilc_projected_potential;
    output_multislice_1 = ilc_projected_potential(system_conf, input_multem);
    toc;
    
    system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 1; % eP_Float = 1, eP_double = 2
    tic;
    clear ilc_projected_potential;
    output_multislice_2 = ilc_projected_potential(system_conf, input_multem);
    toc;
    mean(abs(output_multislice_1.V(:)-output_multislice_2.V(:)))
    
    figure(1);
    subplot(1, 2, 1);
    imagesc(output_multislice_1.V);
    colormap gray;
    axis image;
    subplot(1, 2, 2);
    imagesc(output_multislice_2.V);
    colormap gray;
    axis image;
    disp([min(output_multislice_1.V(:)), min(output_multislice_2.V(:))])
    disp([max(output_multislice_1.V(:)), max(output_multislice_2.V(:))])
    pause(0.10);

end
