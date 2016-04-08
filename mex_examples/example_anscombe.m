% Anscombe transform is a variance-stabilizing transformation 
% that transforms a random variable with a Poisson distribution 
% into one with an approximately standard Gaussian distribution
% https://en.wikipedia.org/wiki/Anscombe_transform
% 
% J = il_anscombe_forward(I) Direct Anscombe transform
% J = il_anscombe_inverse(I) Inverse Anscombe transform
% 
% Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>

clc; clear all;

Im = double(rgb2gray(imread('saturn.png')));

SNR = 5;
[Im_n, k] = il_add_poisson_noise(Im, SNR);

tic;
Im_d = il_anscombe_forward(Im_n);
Im_d = il_filter_wiener(Im_d, 2);
Im_d = il_filter_median(Im_d, 2);
Im_d = il_anscombe_inverse(Im_d);
toc;

figure(1);
subplot(1, 3, 1);
imagesc(Im);
axis image;
colormap gray;
title('original image');

subplot(1, 3, 2);
imagesc(Im_n);
axis image;
colormap gray;
title('Noise image');

subplot(1, 3, 3);
imagesc(Im_d);
axis image;
colormap gray;
title('denoise image');