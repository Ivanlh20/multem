% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath(['..', filesep, 'matlab_functions'])

ilm_mex('release', 'ilc_transf_exp_factor.cpp', '../src');