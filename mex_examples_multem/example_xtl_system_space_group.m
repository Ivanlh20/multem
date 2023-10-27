clear; clc;

addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

% triclinic or anorthic: a
% monoclinic" m
% orthorhombic: o 
% tetragonal: t
% rhombohedral or trigonal: r
% hexagonal: h
% cubic: c

css = 'triclinic';
% css = 'cubic';
css = 'h';
% css = 'o';

% from xtl system string to xtl system number
csn_c = ilc_xtl_css_2_csn(css);
csn_m = ilm_xtl_css_2_csn(css);

disp('from xtl system string to xtl system number')
disp('    C++    Matlab')
disp([csn_c, csn_m])

% from xtl system string to space group range
[sgr_0c, sgr_ec] = ilc_xtl_css_2_sgr(css);
[sgr_0m, sgr_em] = ilm_xtl_css_2_sgr(css);

disp('from xtl system string to space group range')
disp('     C++        Matlab')
disp([sgr_0c, sgr_ec, sgr_0m, sgr_em])
csn = csn_c;

% from xtl system number to space group range
[sgr_0c, sgr_ec] = ilc_xtl_csn_2_sgr(csn);
[sgr_0m, sgr_em] = ilm_xtl_csn_2_sgr(csn);

disp('from xtl system number to space group range')
disp('     C++        Matlab')
disp([sgr_0c, sgr_ec, sgr_0m, sgr_em])

sgn = 125;
% from space group number to xtl system number
csn_c = ilc_xtl_sgn_2_csn(sgn);
csn_m = ilm_xtl_sgn_2_csn(sgn);

disp('from space group number to xtl system number')
disp('    C++    Matlab')
disp([csn_c, csn_m])