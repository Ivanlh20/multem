% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

xtl_parm.na = 1;
xtl_parm.nb = 1;
xtl_parm.nc = 1;

xtl_parm.a = 4.0780;
xtl_parm.b = 4.0780;
xtl_parm.c = 4.0780;
    
xtl_parm.alpha = 90;
xtl_parm.beta = 90;
xtl_parm.gamma = 90;

xtl_parm.sgn = 1;
xtl_parm.pbc = false;

occ = 1;
region = 0;
charge = 0;
% Au = 79
%Z x y z sigma occupancy region charge
rmsd_3d = 0.085;

xtl_parm.asym_uc = [];
    
xtl_parm.base = [79, 0.0, 0.0, 0.0, rmsd_3d, occ, region, charge;...
                    79, 0.5, 0.5, 0.0, rmsd_3d, occ, region, charge;...
                    79, 0.0, 0.5, 0.5, rmsd_3d, occ, region, charge;...
                    79, 0.5, 0.0, 0.5, rmsd_3d, occ, region, charge];

tic;
if 1
    atoms = ilc_xtl_build(xtl_parm);
else
    atoms = ilm_xtl_build(xtl_parm);   
end
toc;

ilm_show_xtl(1, atoms);
