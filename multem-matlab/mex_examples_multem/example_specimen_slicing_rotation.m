% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = multem_input.parameters;         % Load default values;

input_multem.system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_nthread = 4; 
input_multem.system_conf.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multem.simulation_type = 52;	
input_multem.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 110;                      % phonon dimensions (xyz)
input_multem.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multem.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multem.pn_nconf = 100;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multem.spec_rot_theta = 45;                                 % angle (�)
input_multem.spec_rot_u0 = [1 0 0];                               % unitary vector			
input_multem.spec_rot_center_type = 1;                         % 1: geometric center, 2: User define		
input_multem.spec_rot_center_p = [0 0 0];                               % rotation point

na = 8; nb = 8; nc = 8; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.spec_lx = 100;
input_multem.spec_ly = 100;
input_multem.spec_lz = 100;

[input_multem.spec_atoms] = ilm_center_spec(input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, input_multem.spec_lz);

% get spec slicing
[atoms, Slice] = ilc_spec_slicing(input_multem.toStruct);

ilm_show_crystal(1, atoms);

[natoms,~] = size(atoms); 
[nslice, ~] = size(Slice);

for i = 1:nslice
    figure(1); clf;
    i1 = Slice(i, 5); i2 = Slice(i, 6); ii = i1:1:i2;
    plot3(atoms(:, 2), atoms(:, 3), atoms(:, 4), '.k', atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'or');
    set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
    title('Atomic positions');
    ylabel('y','FontSize',14);
    xlabel('x','FontSize',12);
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
        
        xy = [xy; [mean(atoms(ii, 2)), mean(atoms(ii, 3))]];
    end 
end

save('xy_projected.mat', 'xy');
figure(2);
plot(xy(:, 1), xy(:, 2), '*r');