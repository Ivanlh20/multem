% output_multislice = il_MULTEM(system_conf, input_multislice) perform TEM simulation
% 
% Exit wave real space (EWRS) simulation
% 
% All parameters of the input_multislice structure are explained in multem_default_values()
% 
% Copyright 2017 Ivan Lobato <Ivanlh20@gmail.com>

clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;

% eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
% eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
input_multislice.simulation_type = 52;	
input_multislice.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_dim = 110;                      % phonon dimensions (xyz)
input_multislice.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.pn_single_conf = 0;                % 1: true, 0:false (extract single configuration)
input_multislice.pn_nconf = 100;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.spec_rot_theta = 45;                                 % angle (º)
input_multislice.spec_rot_u0 = [1 0 0];                               % unitary vector			
input_multislice.spec_rot_center_type = 1;                         % 1: geometric center, 2: User define		
input_multislice.spec_rot_center_p = [0 0 0];                               % rotation point

na = 8; nb = 8; nc = 8; ncu = 2; rms3d = 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.spec_lx = 100;
input_multislice.spec_ly = 100;
input_multislice.spec_lz = 100;

[input_multislice.spec_atoms] = center_spec(input_multislice.spec_atoms, input_multislice.spec_lx, input_multislice.spec_ly, input_multislice.spec_lz);

% get spec slicing
[atoms, Slice] = il_spec_slicing(input_multislice);

show_crystal(1, atoms);

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