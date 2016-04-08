% [J, k] = il_add_poisson_noise(I, SNR) Add Poisson noise to a signal
%     I is the input signal
%     SNR is the required signal to noise ratio define as
%         'SNR = standard deviation of I / standard deviation of the noise'
%         'standard deviation of the noise = sqrt(variance(J-I))
%     J is the output image with the required SNR
%     k is the scaling factor of I to produce J
% 
% Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>

clc; clear all;

Im = double(rgb2gray(imread('saturn.png')));

SNR = 3;
tic;
[Im_n, k] = il_add_poisson_noise(Im, SNR);
toc;

% calculate output SNR
Im_input = k*Im;
Im_noise = Im_n-Im_input;
SNR_o = std(Im_input(:))/std(Im_noise(:));

disp('output signal to noise ratio')
disp(SNR_o)

figure(1);
subplot(1, 2, 1);
imagesc(Im);
axis image;
colormap gray;
title('original image');

subplot(1, 2, 2);
imagesc(Im_n);
axis image;
colormap gray;
title('noise image');
