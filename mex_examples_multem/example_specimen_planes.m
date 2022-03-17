% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

input_multem = multem_input.parameters;         % Load default values;

input_multem.pn_model = 3;                  % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
input_multem.interaction_model = 1;             % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
input_multem.potential_slicing = 1;             % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
input_multem.pn_dim = 110; 
input_multem.pn_seed = 300183; 
input_multem.pn_nconf = 1;

input_multem.spec_rot_theta = 0;                      % final angle
input_multem.spec_rot_u0 = [1 0 0]; 					% unitary vector			
input_multem.spec_rot_center_type = 1; 			% 1: geometric center, 2: User define		
input_multem.spec_rot_center_p = [0 0 0];					% rotation point

na = 6; nb = 6; nc = 10; ncu = 2; rmsd_3d = 0.085;

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz...
, a, b, c, input_multem.spec_dz] = Au110_xtl(na, nb, nc, ncu, rmsd_3d);

[input_multem.spec_atoms, input_multem.spec_lx...
, input_multem.spec_ly, input_multem.spec_lz] = graphene(10, 2.46, 0.085);

ilm_write_ap_pdb('graphene.pdb', input_multem.spec_atoms);

% occ = 1;
% region = 0;
% charge = 0;
% [input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, lz] = graphene(1, 1.42, sqrt(0.5/(8*pi^2)));
% input_multem.spec_dz = 0.5;

input_multem.spec_dz = 2;

disp([min(input_multem.spec_atoms(:, 4)), max(input_multem.spec_atoms(:,4))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
input_multem.spec_atoms = ilc_add_amorp_lay(input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
toc;

if(lay_pos==1)
    input_multem.spec_amorp(1).z_0 = z_min-lz;          % Starting z position of the amorphous layer (�)
    input_multem.spec_amorp(1).z_e = z_min;             % Ending z position of the amorphous layer (�)
else
    input_multem.spec_amorp(1).z_0 = z_max;             % Starting z position of the amorphous layer (�)
    input_multem.spec_amorp(1).z_e = z_max+lz;          % Ending z position of the amorphous layer (�)
end
input_multem.spec_amorp(1).dz = 2.0;                    % slice thick of the amorphous layer (�)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lz = 10;
Z = 6;
rms_3d = 0.09;
d_min = 1.4;
seed = 1983;
rho = 2.2;
lay_pos = 1; %1: top, 2: bottom

z_min = min(input_multem.spec_atoms(:, 4));
z_max = max(input_multem.spec_atoms(:, 4));

tic;
input_multem.spec_atoms = ilc_add_amorp_lay(input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, lz, d_min, Z, rms_3d, rho, lay_pos, seed);
toc;

if(lay_pos==1)

    input_multem.spec_amorp(2).z_0 = z_min-lz;          % Starting z position of the amorphous layer (�)
    input_multem.spec_amorp(2).z_e = z_min;             % Ending z position of the amorphous layer (�)
else
    input_multem.spec_amorp(2).z_0 = z_max;             % Starting z position of the amorphous layer (�)
    input_multem.spec_amorp(2).z_e = z_max+lz;          % Ending z position of the amorphous layer (�)
end
input_multem.spec_amorp(2).dz = 2.0;                    % slice thick of the amorphous layer (�)

ilm_show_crystal(1, input_multem.spec_atoms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
tic;
[z_planes] = ilc_spec_planes(input_multem.toStruct);
toc;
diff(z_planes)
[nplanes, ~] = size(z_planes);
disp(['Number of slices = ', num2str(nplanes)])

figure(1); clf;
plot(input_multem.spec_atoms(:, 3), input_multem.spec_atoms(:, 4), 'ok');
set(gca,'ydir','reverse');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Atomic positions');
ylabel('y','FontSize',14);
xlabel('x','FontSize',12);
axis equal;

for i = 1:nplanes
    hold on;
    plot([-2 input_multem.spec_lx], [z_planes(i) z_planes(i)], '-r');    
    axis equal;
end
axis([-2, 18, min(input_multem.spec_atoms(:, 4))-5, max(input_multem.spec_atoms(:, 4))+5]);

% nbins = floor((max(input_multem.spec_atoms(:, 4))-min(input_multem.spec_atoms(:, 4)))/0.10001);
% tic;
% [x, y] = ilc_hist(input_multem.spec_atoms(:, 4), nbins-1);
% toc;
% 
% figure(2); clf;
% plot(x, y, '-+r');
% hold on;
% ii = find(y<0.5);
% plot(x(ii), y(ii), '.b');