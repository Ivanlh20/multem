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

na = 4; nb = 4; nc = 4; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 2048; 
input_multislice.ny = 2048;

clear get_specimen_slicing;
[atoms, Slice] = get_specimen_slicing(input_multislice);

[natoms,~] = size(atoms); [nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multislice.islice = islice;
    
    input_multislice.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    input_multislice.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear get_projected_potential;
    ouput_multislice_1 = get_projected_potential(input_multislice);
    toc;
    
    input_multislice.device = 1;                        % eD_CPU = 1, eD_GPU = 2
    input_multislice.precision = 2;                     % eP_Float = 1, eP_double = 2
    tic;
    clear get_projected_potential;
    ouput_multislice_2 = get_projected_potential(input_multislice);
    toc;
    sum(abs(ouput_multislice_1.V(:)-ouput_multislice_2.V(:))/(input_multislice.nx*input_multislice.ny))
    
    figure(1);   
    imagesc(ouput_multislice_1.V);
    colormap gray;
    axis image;  
    figure(2);   
    imagesc(ouput_multislice_2.V);
    colormap gray;
    axis image; 
    
    pause(0.10);
end;
