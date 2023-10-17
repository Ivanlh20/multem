% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath([ fileparts(pwd), filesep, 'matlab_functions'])

ilm_mex('release', 'ilc_iehwgd_2_sigma.cpp', '../src');
