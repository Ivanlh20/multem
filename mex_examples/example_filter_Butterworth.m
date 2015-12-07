clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Radius = nx/4; n = 16; shift = 0;
tic;
f_1d = il_filter_Butterworth(ny, 1, dy, dx, Radius, n, shift);
toc;

tic;
f_2d = il_filter_Butterworth(ny, nx, dy, dx, Radius, n, shift);
toc;

figure(1);
subplot(1, 2, 1);
plot(f_1d);
subplot(1, 2, 2);
imagesc(f_2d);
colormap hot;
axis image;