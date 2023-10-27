% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

E_0 = 300;
bs = [50, 35]; % simulations box size
theta = [80, 100, 120, 200, 215, 250]; % mrad

% minimun number of pixels to cover the detector
np = ilc_det_min_spl(E_0, bs, theta);
disp(np)