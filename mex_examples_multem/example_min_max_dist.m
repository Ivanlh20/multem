% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

[x, y, z] = meshgrid(0:2:50);
r_3d = [x(:), y(:), z(:)];
r_3d = r_3d + 0.1*rand(size(r_3d));
r_3d = double(r_3d);

tic;
[r, rdf] = ilc_rdf(r_3d, 5, 400);
toc;

r_min = ilc_min_dist(r_3d);
r_max = ilc_max_dist(r_3d);

figure(1); clf;
plot(r, rdf, '-k');
hold on;
plot(r_min, 0.5*max(rdf), 'or')
hold on;
plot(r_max, 0.5*max(rdf), 'or')