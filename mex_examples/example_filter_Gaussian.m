clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Sigma = 128; shift = 0;

tic;
f_1d = get_filter_Gaussian(ny, 1, dy, dx, Sigma, shift);
toc;

tic;
f_2d = get_filter_Gaussian(ny, nx, dy, dx, Sigma, shift);
toc;

figure(1);
subplot(1, 2, 1);
plot(f_1d);
subplot(1, 2, 2);
imagesc(f_2d);
colormap hot;
axis image;