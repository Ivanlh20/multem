clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Sigma = 128; shift = 0;

tic;
f_1d = il_func_gaussian(nx, 1, dx, dy, Sigma, shift);
toc;

tic;
f_2d = il_func_gaussian(nx, ny, dx, dy, Sigma, shift);
toc;

tic;
f_2d_r = il_func_gaussian_by_row(nx, ny, dx, dy, Sigma, shift);
toc;

figure(1);
subplot(1, 3, 1);
plot(f_1d, '-r');
xlim([0 nx]);
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1 1 1]);
title('1d Gaussian function')

subplot(1, 3, 2);
imagesc(f_2d);
colormap hot;
axis image;
title('2d Gaussian function')

subplot(1, 3, 3);
imagesc(f_2d_r);
colormap hot;
axis image;
title('2d Gaussian function by row')