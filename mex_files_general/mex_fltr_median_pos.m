% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
path_add = ['..', filesep, 'matlab_functions'];

ilm_mex('release', 'ilc_fltr_median_pos.cpp', '../src');