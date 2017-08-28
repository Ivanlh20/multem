clear all; clc;

input_multislice = multem_default_values();         % Load default values;

system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
system_conf.cpu_nthread = 4; 
system_conf.gpu_device = 0;

input_multislice.pn_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 2;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.potential_type = 6;                % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

input_multislice.pn_dim = 111;                      % phonon dimensions
input_multislice.pn_seed = 300183;                  % Random seed(frozen phonon)
input_multislice.pn_single_conf = 1;                % 1: true, 0:false
input_multislice.pn_nconf = 1;                      % true: phonon configuration, false: number of frozen phonon configurations

na = 4; nb = 4; nc = 4; ncu = 2; rms3d = 0.085;

[input_multislice.spec_atoms, input_multislice.spec_lx...
, input_multislice.spec_ly, input_multislice.spec_lz...
, a, b, c, input_multislice.spec_dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.nx = 2048; 
input_multislice.ny = 2048;

clear il_spec_slicing;
[atoms, Slice] = il_spec_slicing(input_multislice);

[natoms,~] = size(atoms); [nslice, ~] = size(Slice);
for islice = 1:nslice
    input_multislice.islice = islice;
    
    system_conf.device = 1;                        % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
    tic;
    clear il_projected_potential;
    ouput_multislice_1 = il_projected_potential(system_conf, input_multislice);
    toc;
    
    system_conf.device = 2;                        % eD_CPU = 1, eD_GPU = 2
    system_conf.precision = 1;                     % eP_Float = 1, eP_double = 2
    tic;
    clear il_projected_potential;
    ouput_multislice_2 = il_projected_potential(system_conf, input_multislice);
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
