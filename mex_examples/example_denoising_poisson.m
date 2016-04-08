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

% add poisson noise
SNR = 5; % signal to noise ratio
tic;
[Im_n, k] = il_add_poisson_noise(Im, SNR);
Im = k*Im;
toc;

figure(1);
subplot(1, 4, 1);
imagesc(Im);
axis image;
colormap gray;
title('original image');

subplot(1, 4, 2);
imagesc(Im_n);
axis image;
colormap gray;
title('Noise image');

tic;
Im_pd = il_denoising_poisson(Im_n, 2, 2);
toc;
subplot(1, 4, 3);
imagesc(Im_pd);
axis image;
colormap gray;
title('denoise image');

tic;
Lr = il_denoising_poisson_by_row(Im_n, 2, 2);
toc;
subplot(1, 4, 4);
imagesc(Lr);
axis image;
colormap gray;
title('denoise image by row');

figure(2);
x = 1:size(Im, 2);
plot(x, Im_n(700, :), '-k', x, Im(700, :), '-b', x, Im_pd(700, :), '-r');
