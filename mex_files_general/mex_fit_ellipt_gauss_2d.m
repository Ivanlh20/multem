% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath(['..', filesep, 'matlab_functions'])

ilm_mex('release', 'ilc_fit_ellipt_gauss_2d.cpp', '../src');