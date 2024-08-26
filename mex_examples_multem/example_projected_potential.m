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

input_multem.pn_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 2;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 111;                      % phonon dimensions
input_multem.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multem.pn_single_conf = 1;                % 1: true, 0:false
input_multem.pn_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

na = 4; nb = 4; nc = 4; ncu = 2; rmsd_3d = 0.25;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);
input_multem.spec_dz = 10000;
input_multem.nx = 2048; 
input_multem.ny = 2048;

clear ilc_spec_slicing;
[atoms, Slice] = ilc_spec_slicing(input_multem.toStruct);

[natoms,~] = size(atoms); [nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multem.islice = islice;
    
    input_multem.system_conf.device = 1;                        % eD_CPU = 1, eD_GPU = 2
    input_multem.system_conf.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear ilc_projected_potential;
    ouput_multislice_1 = input_multem.ilc_projected_potential;
    toc;
    
    input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    input_multem.system_conf.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear ilc_projected_potential;
    ouput_multislice_2 = input_multem.ilc_projected_potential;
    toc;
    mean(abs(ouput_multislice_1.V(:)-ouput_multislice_2.V(:)))
    
    figure(1); 
    subplot(1, 2, 1);
    imagesc(ouput_multislice_1.V);
    colormap gray;
    axis image;  
    subplot(1, 2, 2);  
    imagesc(ouput_multislice_2.V);
    colormap gray;
    axis image; 
    disp([min(ouput_multislice_1.V(:)), min(ouput_multislice_2.V(:))])
    disp([max(ouput_multislice_1.V(:)), max(ouput_multislice_2.V(:))])
    pause(0.10);

end
