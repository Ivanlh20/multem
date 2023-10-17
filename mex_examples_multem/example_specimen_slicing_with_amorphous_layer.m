clear;clc;

input_multem = ilm_dflt_input_multem(); % Load default values;

input_multem.atomic_vib_mod = 1; % eavm_still_atom = 1, eavm_absorptive_pot = 2, eavm_frozen_phonon = 3, eavm_user_def = 4
input_multem.elec_spec_interact_mod = 1; % eesim_multislice = 1, eesim_phase_object = 2, eesim_weak_phase_object = 3
input_multem.spec_slic(1).typ = 1; % esst_plns_proj = 1, esst_dz_proj = 2, esst_plns_sub = 3, esst_dz_sub = 4, esst_user_def = 5, esst_auto = 6
input_multem.atomic_vib_dim = [true, true, false];
input_multem.atomic_vib_seed = 300183;
input_multem.atomic_vib_nconf = 1;

input_multem.spec_rot_theta = 0; % final angle
input_multem.spec_rot_u_0 = [1 0 0]; % unitary vector			
input_multem.spec_rot_ctr_typ = 1; % 1: geometric center, 2: User define		
input_multem.spec_rot_ctr_p = [0 0 0]; % rotation point

na = 6;nb = 6;nc = 10;ncu = 4;rmsd_3d = 0.15;

[input_multem.spec_atoms, input_multem.spec_bs_x...
, input_multem.spec_bs_y, input_multem.spec_bs_z...
, a, b, c, input_multem.spec_slic(1).sli_thick] = Au001_xtl(na, nb, nc, ncu, rmsd_3d);

input_multem.spec_slic(1).sli_thick=a/2;

disp([min(input_multem.spec_atoms(:, 4)), max(input_multem.spec_atoms(:, 4))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lz = 20;
Z = 6;
rms_3d = 0.09;
d_min = 1.4;
seed = 1983;
rho = 2.2;
lay_pos = 2; %1: top, 2: bottom

z_min = min(input_multem.spec_atoms(:, 4));
z_max = max(input_multem.spec_atoms(:, 4));
tic;
input_multem.spec_atoms = ilc_amorp_lay_add(input_multem.spec_atoms, input_multem.spec_bs_x, input_multem.spec_bs_y, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
toc;

if(lay_pos==1)
    input_multem.spec_slic(1).z_0 = z_min-lz; % Starting z position of the amorphous layer (Å)
    input_multem.spec_slic(1).z_e = z_min; % Ending z position of the amorphous layer (Å)
else
    input_multem.spec_slic(1).z_0 = z_max; % Starting z position of the amorphous layer (Å)
    input_multem.spec_slic(1).z_e = z_max+lz; % Ending z position of the amorphous layer (Å)
end
input_multem.spec_slic(1).sli_thick = 2.0; % slice thick of the amorphous layer (Å)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lz = 10;
% Z = 6;
% rms_3d = 0.09;
% d_min = 1.4;
% seed = 1983;
% rho = 2.2;
% lay_pos = 1; %1: top, 2: bottom
% 
% z_min = min(input_multem.spec_atoms(:, 4));
% z_max = max(input_multem.spec_atoms(:, 4));
% 
% tic;
% input_multem.spec_atoms = ilc_amorp_lay_add(input_multem.spec_atoms, input_multem.spec_bs_x, input_multem.spec_bs_y, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
% toc;
% 
% if(lay_pos==1)
%  input_multem.spec_slic(2).z_0 = z_min-lz; % Starting z position of the amorphous layer (Å)
%  input_multem.spec_slic(2).z_e = z_min; % Ending z position of the amorphous layer (Å)
% else
%  input_multem.spec_slic(2).z_0 = z_max; % Starting z position of the amorphous layer (Å)
%  input_multem.spec_slic(2).z_e = z_max+lz; % Ending z position of the amorphous layer (Å)
% end
% input_multem.spec_slic(2).dz = 2.0; % slice thick of the amorphous layer (Å)

% ilm_show_xtl(1, input_multem.spec_atoms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
tic;
[atoms, Slice] = ilc_spec_slicing(input_multem);
toc;
disp([min(atoms(:, 4)), max(atoms(:, 4))])

[nslice, ~] = size(Slice);
disp(['Number of slices = ', num2str(nslice)])

figure(1); clf;
plot(atoms(:, 3), atoms(:, 4), 'ok');
set(gca, 'ydir', 'reverse');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
title('Atomic positions');
ylabel('y', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
axis equal;

for i = 1:nslice
    hold on;
    plot([-2 input_multem.spec_bs_x], [Slice(i, 1) Slice(i, 1)], '-r', [-2 input_multem.spec_bs_x], [Slice(i, 2) Slice(i, 2)], '-r');
    axis equal;

end
axis([-2, 18, min(input_multem.spec_atoms(:, 4))-5, max(input_multem.spec_atoms(:, 4))+5]);

tic;
[z_planes] = ilc_spec_planes(input_multem);
toc;
nplanes = length(z_planes);
for i = 1:nplanes
    hold on;
    plot([-2 input_multem.spec_bs_x], [z_planes(i) z_planes(i)], '-b');
    axis equal;
end
diff(z_planes)

nbins = floor((max(atoms(:, 4))-min(atoms(:, 4)))/0.10001);
tic;
[x, y] = ilc_hist(input_multem.spec_atoms(:, 4), nbins-1);
toc;

% figure(2);clf;
% plot(x, y, '-+r');
% hold on;
% ii = find(y<0.5);
% plot(x(ii), y(ii), '.b');