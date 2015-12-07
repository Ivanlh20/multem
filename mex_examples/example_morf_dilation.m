clear all; clc;

Im = double(rgb2gray(imread('peppers.png')));

figure(1) 
subplot(1, 3, 1);
imagesc(Im);
axis image;
colormap gray;

subplot(1, 3, 2);
level = il_otsu_threshold(Im, 256);
BW = il_binary(Im, level);
imagesc(BW);
axis image;
colormap gray;

subplot(1, 3, 3);
tic;
BW_d = il_morf_dilation(BW, 1);
toc;
imagesc(BW_d-BW);
axis image;
colormap gray;

