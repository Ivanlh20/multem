clc; clear all;
addpath( '../matlab_functions')
  
ilm_mex('release', 'ilc_projected_potential.cu', '../src');