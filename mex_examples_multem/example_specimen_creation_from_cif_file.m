% Specimen creation
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear all;clc;

ncu = 2;
rmsd_3d = 0.085;
fn = 'SrTiO3_mp-4651_conventional_standard.cif';

[a, b, c] = ilm_read_lat_parm_cif(fn);
na = 2; 
nb = 2; 
nc = 2;
rmsd_3d_0 = 0.085;
pbc = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[atoms_0, lx, ly, lz] = ilm_read_ap_cif(fn, rmsd_3d_0, pbc, na, nb, nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbc = true;
[atoms_1, lx, ly, lz] = ilm_read_ap_cif(fn, rmsd_3d_0, pbc, na, nb, nc);

figure(1); clf;
subplot(1, 2, 1);
ilm_show_crystal(0, atoms_0, false);
title('pbc=false');
subplot(1, 2, 2);
ilm_show_crystal(0, atoms_1, false);
title('pbc=true');
