clc; clear all;
addpath( '../matlab_functions')

ilm_mex('release', 'ilc_fwhm_2_sigma.cpp', '../src');
