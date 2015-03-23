clear all;
clc;
m = 2e+8;
z = rand(m, 1);
% Mex file single core
tic; get_QuickSort_CPU(z); toc;
% Matlab function mul-core
tic; sort(z); toc;