clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; Radius = nx/4; n = 16; shift = 0;

tic;
f_1d = il_func_butterworth(nx, 1, dx, dy, Radius, n, shift);
toc;

tic;
f_2d = il_func_butterworth(nx, ny, dx, dy, Radius, n, shift);
toc;

tic;
f_2d_r = il_func_butterworth_by_row(nx, ny, dx, dy, Radius, n, shift);
toc;

figure(1);
subplot(1, 3, 1);
plot(f_1d, '-r');
xlim([0 nx]);
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1 1 1]);
title('1d Butterworth function')

subplot(1, 3, 2);
imagesc(f_2d);
colormap hot;
axis image;
title('2d Butterworth function')

subplot(1, 3, 3);
imagesc(f_2d_r);
colormap hot;
axis image;
title('2d Butterworth function by row')