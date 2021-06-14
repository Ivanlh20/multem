% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

xtl_parm.a = 5;
xtl_parm.b = 5;
xtl_parm.c = 5;
xtl_parm.alpha = 90;
xtl_parm.beta = 90;
xtl_parm.gamma = 90;
xtl_parm.na = 1;
xtl_parm.nb = 1;
xtl_parm.nc = 1;
xtl_parm.sgn = 120;
xtl_parm.pbc = false;
xtl_parm.asym_uc = [6, 0, 0, 0; 6, 1/4, 1/4, 1/4];
if 1
    xtl_parm.base = ilc_xtl_build_base(xtl_parm.asym_uc, xtl_parm.sgn);
else
    xtl_parm.base = ilm_xtl_build_base(xtl_parm); 
end

base = xtl_parm.base;

ilm_show_xtl(1, base);
view([0 0 1])