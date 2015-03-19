clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Radius = 128; n = 4; shift = 1;
f = getButterworth_Filter_2D(ny, nx, dy, dx, Radius, n, shift);
figure(1);
imagesc(f);
colormap gray;
axis image;