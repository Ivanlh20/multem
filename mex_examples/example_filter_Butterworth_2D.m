clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Radius = nx/4; n = 16; shift = 0;
tic;
f = get_filter_Butterworth_2D(ny, nx, dy, dx, Radius, n, shift);
toc;

figure(1);
imagesc(f);
colormap hot;
axis image;