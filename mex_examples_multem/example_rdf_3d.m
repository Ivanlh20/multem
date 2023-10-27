% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'mex_bin'])

[x, y, z] = meshgrid(0:2:50);
r_3d = [x(:), y(:), z(:)];
r_3d = r_3d + 0.1*rand(size(r_3d));
r_3d = double(r_3d);

tic;
[r, rdf] = ilc_rdf(r_3d, 8, 200);
toc;

figure(1); clf;
plot(r, rdf, '-+r');