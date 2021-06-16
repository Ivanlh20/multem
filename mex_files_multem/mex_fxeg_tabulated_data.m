% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath(['..', filesep, 'matlab_functions'])

ilm_mex('debug', 'ilc_fxeg_tabulated_data.cpp', '../src');