% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clc; clear;
addpath([ fileparts(pwd), filesep, 'matlab_functions'])
  
ilm_mex('release', 'ilc_elec_interact_parm_kva.cpp', '../src');