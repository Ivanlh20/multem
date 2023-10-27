% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath(['..', filesep, 'matlab_functions'])

ilm_mex('release', 'ilc_chi2_pcf_2d.cu', '../src');