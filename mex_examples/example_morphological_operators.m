clear all; clc;

Im = double(rgb2gray(imread('peppers.png')));

figure(1) 
subplot(2, 2, 1);
imagesc(Im);
axis image;
colormap gray;
title('Image');

subplot(2, 2, 2);
level = il_otsu_threshold(Im, 256);
tic;
BW = il_binarization(Im, level);
toc;
imagesc(BW);
axis image;
colormap gray;
title('Otsu threshold');

BW = Im;

subplot(2, 2, 3);
tic;
BW_d = il_morp_erode(BW, 1);
toc;
imagesc(BW_d);
axis image;
colormap gray;
title('erode');

subplot(2, 2, 4);
tic;
BW_d = il_morp_dilate(BW, 1);
toc;
imagesc(BW_d);
axis image;
colormap gray;
title('dilate');

