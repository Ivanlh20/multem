clc; clear;
addpath(['..', filesep, 'matlab_functions'])

ilm_mex('release', 'ilc_xtl_build_base.cpp', '../src');