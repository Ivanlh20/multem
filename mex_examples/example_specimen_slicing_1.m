clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.phonon_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.fp_dim = 111; 
input_multislice.fp_seed = 300183; 
input_multislice.fp_nconf = 3;

input_multislice.tm_active = 0;						% 1: true, 0:false
input_multislice.tm_theta = 5; 					% final angle
input_multislice.tm_u0 = [1 0 0]; 					% unitary vector			
input_multislice.tm_rot_point_type = 1; 			% 1: geometric center, 2: User define		
input_multislice.tm_p0 = [0 0 0];					% rotation point

na = 4; nb = 4; nc = 20; ncu = 2; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

input_multislice.dz = 5.0;

% get specimen slicing
tic;
[atoms, Slice] = il_specimen_slicing(input_multislice);
toc;
[natoms,~] = size(atoms); [nslice, ~] = size(Slice);

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
end;

[size(input_multislice.atoms, 1), natoms, nslice]
[input_multislice.lx, input_multislice.ly, input_multislice.lz]