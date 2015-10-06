clear all; clc;

input_multislice = multem_default_values();         % Load default values;

input_multislice.phonon_model = 1;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 3;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.fp_dim = 110; 
input_multislice.fp_seed = 300183; 
input_multislice.fp_nconf = 1;

input_multislice.tm_active = 0;						% 1: true, 0:false
input_multislice.tm_nrot = 10; 						% number of rotations
input_multislice.tm_irot = 1;						% specific rotation configuration
input_multislice.tm_theta_0 = -10; 					% initial angle
input_multislice.tm_theta_e = +10; 					% final angle
input_multislice.tm_u0 = [0 1 1]; 					% unitary vector			
input_multislice.tm_rot_point_type = 1; 			% 1: geometric center, 2: User define		
input_multislice.tm_p0 = [0 0 0];					% rotation point

input_multislice.lx = 10;
input_multislice.ly = 10;
input_multislice.lz = 10;
input_multislice.dz = 0.5;

input_multislice.atoms = [29, 2, 2, 0.0, 0.8, 1.0; 29, 6, 2, 0.0, 0.8, 1.0];
[input_multislice.atoms, input_multislice.lx, input_multislice.ly, lz] = graphene(1, sqrt(0.5/(8*pi^2)));
input_multislice.dz = 0.5;

% get specimen slicing
tic;
input_multislice.phonon_model = 1;
[atoms0, Slice0] = get_specimen_slicing(input_multislice);
toc;

[nslice0, ~] = size(Slice0);

tic;
input_multislice.phonon_model = 3;
[atoms, Slice] = get_specimen_slicing(input_multislice);
toc;

[nslice, ~] = size(Slice);

figure(1); clf;
plot(atoms(:, 2), atoms(:, 4), '*k');   
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Atomic positions');
ylabel('y','FontSize',14);
xlabel('x','FontSize',12);
axis equal;
axis([-2 input_multislice.lx+2 -5 input_multislice.lz + 5]);


for i = 1:nslice
    hold on;
    plot([-2 18], [Slice(i, 1) Slice(i, 1)], '-b', [-2 18], [Slice(i, 2) Slice(i, 2)], '-b');    
    axis equal;
    axis([-2 input_multislice.lx+2 -5 input_multislice.lz + 5]);
end;

for i = 1:nslice0
    hold on;
    plot([-2 input_multislice.lx+2], [Slice0(i, 1) Slice0(i, 1)], '-r', [-2 input_multislice.lx+2], [Slice0(i, 2) Slice0(i, 2)], '-r');    
    axis equal;
    axis([-2 input_multislice.lx+2 -5 input_multislice.lz + 5]);
end;