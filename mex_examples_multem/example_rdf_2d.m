% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'mex_bin'])

[x, y] = meshgrid(0:2:50);
r_2d = [x(:), y(:)];
r_2d = r_2d + 0.1*rand(size(r_2d));
r_2d = double(r_2d);

tic;
[r, rdf] = ilc_rdf(r_2d, 8, 200);
toc;

figure(1); clf;
plot(r, rdf, '-+r');