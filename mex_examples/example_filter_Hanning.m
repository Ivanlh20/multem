clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; k = 0.5; shift = 0;
tic;
f_1d = il_filter_Hanning(ny, 1, dy, dx, k, shift);
toc;

tic;
f_2d = il_filter_Hanning(ny, nx, dy, dx, k, shift);
toc;

figure(1);
subplot(1, 2, 1);
plot(f_1d);
subplot(1, 2, 2);
imagesc(f_2d);
colormap hot;
axis image;