% output_multem = ilc_multem(system_conf, input_multem) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multem structure are explained in ilm_dflt_input_multem()
% 
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

system_conf.precision = 1; % eP_Float = 1, eP_double = 2
system_conf.device = 2; % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_n_thread = 4;
system_conf.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_STEM_EELS=61, eTEMST_ISTEM_EELS=62, eTEMST_EFTEMFS=71, eTEMST_EFTEMRS=72, eTEMST_ProbeFS=81, eTEMST_ProbeRS=82, eTEMST_PPFS=91, eTEMST_PPRS=92, eTEMST_TFFS=101, eTEMST_TFRS=102
input_multem.simulation_type = 52;
input_multem.pn_model = 1; % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1; % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1; % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6; % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 110; % phonon dimensions (xyz)
input_multem.pn_seed = 300183; % Random seed(frozen phonon)
input_multem.pn_single_conf = 0; % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 100; % true: phonon configuration, false: number of frozen phonon configurations

input_multem.spec_rot_theta = 45; % angle (º)
input_multem.spec_rot_u0 = [1 0 0]; % unitary vector			
input_multem.spec_rot_center_type = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_center_p = [0 0 0]; % rotation point

na = 8;nb = 8;nc = 8;ncu = 2;rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.spec_bs_x = 100;
input_multem.spec_bs_y = 100;
input_multem.spec_bs_z = 100;

[input_multem.spec_atoms] = ilm_spec_recenter(input_multem.spec_atoms, input_multem.spec_bs_x, input_multem.spec_bs_y, input_multem.spec_bs_z);

% get spec slicing
[atoms, Slice] = ilc_spec_slicing(input_multem);

ilm_show_xtl(1, atoms);

[natoms, ~] = size(atoms);
[nslice, ~] = size(Slice);

for i = 1:nslice
    figure(1); clf;
    i1 = Slice(i, 5);i2 = Slice(i, 6);ii = i1:1:i2;
    plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), '.k', atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'or');
    set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
    title('Atomic positions');
    ylabel('y', 'FontSize', 14);
    xlabel('x', 'FontSize', 12);
    axis equal;
    i2-i1+1
    view([1 0 0]);
    pause(0.1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
natoms = size(atoms, 1);
bb = zeros(natoms, 1);
d = 0.1;
ic = 0;
xy = [];
for ia=1:natoms
    if(bb(ia)<0.1)
        x = atoms(ia, 2);
        y = atoms(ia, 3);
        ii = find(sqrt((atoms(:, 2)-x).^2+(atoms(:, 3)-y).^2)<d);
        bb(ii) = 1;
        
        xy = [xy;[mean(atoms(ii, 2)), mean(atoms(ii, 3))]];
    end 
end

save('xy_projected.mat', 'xy');
figure(2);
plot(xy(:, 1), xy(:, 2), '*r');