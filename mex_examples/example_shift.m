clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Sigma = 1; shift = 0;

tic;
f = il_filter_Gaussian(ny, nx, dy, dx, Sigma, shift);
toc;

x_shift = -256;
y_shift = 256;

tic;
f_shift = il_shift(f, dy, dx, x_shift, y_shift);
toc;
figure(1);
subplot(1, 2, 1);
imagesc(f);
colormap hot;
axis image;
subplot(1, 2, 2);
imagesc(f_shift);
colormap hot;
axis image;
