% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath([ fileparts(pwd), filesep, 'matlab_functions'])

ilm_mex('debug', 'ilc_amorp_build.cpp', '../src');
