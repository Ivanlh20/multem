clc; clear all;
addpath( '../matlab_functions')

ilm_mex('release', 'ilc_propagate.cu', '../src');