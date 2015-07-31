clear all; clc;

input_multislice.phonon_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multislice.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multislice.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multislice.fp_dim = 111; 
input_multislice.fp_seed = 300183; 
input_multislice.fp_iconf = 1;

na = 4; nb = 4; nc = 10; ncu = 4; rms3d = 0.085;

[input_multislice.atoms, input_multislice.lx...
, input_multislice.ly, input_multislice.lz...
, a, b, c, input_multislice.dz] = Au001Crystal(na, nb, nc, ncu, rms3d);

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
axis([-2 18 -5 input_multislice.lz + 5]);


for i = 1:nslice
    hold on;
    plot([-2 18], [Slice(i, 1) Slice(i, 1)], '-b', [-2 18], [Slice(i, 2) Slice(i, 2)], '-b');    
    axis equal;
    axis([-2 18 -5 input_multislice.lz + 5]);
end;

for i = 1:nslice0
    hold on;
    plot([-2 18], [Slice0(i, 1) Slice0(i, 1)], '-r', [-2 18], [Slice0(i, 2) Slice0(i, 2)], '-r');    
    axis equal;
    axis([-2 18 -5 input_multislice.lz + 5]);
end;