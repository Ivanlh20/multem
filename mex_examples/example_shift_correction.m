clc; clear all;
nx = 1024; 
ny = 1024;
lx = 30;
ly = 30;

dx = lx/nx; dy = ly/ny; Sigma = 1; shift = 0;
SNR = 5; % signal to noise ratio

% create a gaussian function and add noise 
tic;
Im = il_func_gaussian(nx, ny, dx, dy, Sigma, shift);
Im_0 = il_add_poisson_noise(Im, SNR);
toc;

xs_0 = -5;
ys_0 = 7;

% shift image and add noise
tic;
Im_s = il_shift(Im, dx, dy, xs_0, ys_0);
Im_s = il_add_poisson_noise(Im_s, SNR);
toc;

% find drift
k = 0.1;
sigma_cut_off = (nx*dx/2)/2;

% calculate pcf
tic;
pcf = il_pcf(Im_0, Im_s, dx, dy, k, sigma_cut_off);
toc;

% % find shift and corrected
% tic;
% [xs, ys] = il_find_shift_2d(Im_0, Im_s, dx, dy, k, sigma_cut_off);
% Im_s = il_shift(Im_s, dx, dy, -xs, -ys);
% toc;
% [xs, ys]

% correct shift
tic;
[Im_sc, xs, ys] = il_correct_shift_2d(Im_0, Im_s, dx, dy, k, sigma_cut_off);
toc;
disp([xs, ys]);

figure(1);
subplot(2,2, 1);
imagesc(Im_0);
colormap hot;
axis image;
title('input image')

subplot(2, 2, 2);
imagesc(Im_s);
colormap hot;
axis image;
title('shift image')

subplot(2, 2, 3);
imagesc(Im_sc);
colormap hot;
axis image;
title('shift corretion')

subplot(2, 2, 4);
imagesc(pcf);
colormap hot;
axis image;
title('Phase correlation')