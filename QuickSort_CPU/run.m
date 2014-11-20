clear all;
clc;
m = 2e+8;
z = rand(m, 1);
% Mex file
tic;
zs = QuickSort(z);
toc;
% Matlab
tic;
zs = sort(z);
toc;
tic
a = zs;
toc;