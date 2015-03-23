clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Sigma = 128; shift = 1;
f = get_Gaussian_Filter_2D_CPU(ny, nx, dy, dx, Sigma, shift);
figure(1);
imagesc(f);
colormap gray;
axis image;