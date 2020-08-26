clc; clear all;
addpath( '../matlab_functions')

ilm_mex('release', 'ilc_hwhm_2_sigma.cpp', '../src');
