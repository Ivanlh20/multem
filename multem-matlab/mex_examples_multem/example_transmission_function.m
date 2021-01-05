% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = multem_input.parameters;         % Load default values;

input_multem.system_conf.precision = 1;                          % eP_Float = 1, eP_double = 2
input_multem.system_conf.device = 2;                             % eD_CPU = 1, eD_GPU = 2
input_multem.system_conf.cpu_ncores = 1;                         % Number of Cores CPU (It will be used in the future)
input_multem.system_conf.cpu_nthread = 4;                   % Number of CPU threads 
input_multem.system_conf.gpu_device = 0;                    % GPU device (i.e. 0, 1, 2, ... )

input_multem.pn_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multem.pn_dim = 110;                      % phonon dimensions
input_multem.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multem.pn_single_conf = 1;                % 1: true, 0:false
input_multem.pn_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multem.bwl = 0;

input_multem.E_0 = 300;
input_multem.theta = 0.0; 
input_multem.phi = 0.0;

na = 4; nb = 4; nc =2; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

% input_multem.spec_atoms = [79 4.0 4.0 0 rmsd_3d 1.0];
% input_multem.spec_lx = 8.0; 
% input_multem.spec_ly = 8.0; 
% input_multem.spec_dz = 0.5;

input_multem.nx = 2048; 
input_multem.ny = 2048;

clear ilc_spec_slicing;
[atoms, Slice] = ilc_spec_slicing(input_multem.toStruct);

[natoms,~] = size(atoms); [nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multem.islice = islice;
    
    input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    input_multem.system_conf.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
%     clear ilc_transmission_function;
    ouput_multislice_1 = input_multem.ilc_transmission_function;
    toc;
    
    input_multem.system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    input_multem.system_conf.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
%     clear ilc_transmission_function;
    ouput_multislice_2 = input_multem.ilc_transmission_function;
    toc;
    sum(abs(ouput_multislice_1.trans(:)-ouput_multislice_2.trans(:))/(input_multem.nx*input_multem.ny))
    
    figure(1);
    subplot(1, 3, 1);    
    imagesc(real(ouput_multislice_1.trans));
    colormap gray;
    axis image;
    subplot(1, 3, 2);    
    imagesc(imag(ouput_multislice_1.trans));
    colormap gray;
    axis image;   
    subplot(1, 3, 3);    
    imagesc(abs(ouput_multislice_1.trans));
    colormap gray;
    axis image;   
    num2str([islice, min(abs(ouput_multislice_1.trans(:))), max(abs(ouput_multislice_1.trans(:))), sum(abs(ouput_multislice_1.trans(:)))/(input_multem.nx*input_multem.ny)], 10)
    pause(0.10);
end