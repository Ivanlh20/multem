clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.precision = 1;                     % eP_Float = 1, eP_double = 2
input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
input_multislice.cpu_ncores = 1; 
input_multislice.cpu_nthread = 4; 
input_multislice.gpu_device = 0;
input_multislice.gpu_nstream = 8;

input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.fp_dim = 111;                      % phonon dimensions
input_multislice.fp_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.fp_single_conf = 1;                % 1: true, 0:false
input_multislice.fp_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

input_multislice.bwl = 1;

input_multislice.E_0 = 300;
input_multislice.theta = 0.0; 
input_multislice.phi = 0.0;

na = 4; nb = 4; nc = 5; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

% input_multislice.atoms = [79 4.0 4.0 0 rms3d 1.0];
% input_multislice.lx = 8.0; 
% input_multislice.ly = 8.0; 
% input_multislice.dz = 0.5;

input_multislice.nx = 1024; 
input_multislice.ny = 1024;

clear get_specimen_slicing;
[atoms, Slice] = get_specimen_slicing(input_multislice);

[natoms,~] = size(atoms); [nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multislice.islice = islice;
    
    input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
    input_multislice.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear get_transmission_function;
    ouput_multislice_1 = get_transmission_function(input_multislice);
    toc;
    
    input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    input_multislice.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear get_transmission_function;
    ouput_multislice_2 = get_transmission_function(input_multislice);
    toc;
    sum(abs(ouput_multislice_1.trans(:)-ouput_multislice_2.trans(:))/(input_multislice.nx*input_multislice.ny))
    
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
    num2str([islice, min(abs(ouput_multislice_1.trans(:))), max(abs(ouput_multislice_1.trans(:))), sum(abs(ouput_multislice_1.trans(:)))/(input_multislice.nx*input_multislice.ny)], 10)
    pause(0.10);
end;