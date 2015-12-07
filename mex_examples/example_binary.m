clear all; clc;

Im = double(rgb2gray(imread('peppers.png')));

figure(1) 
subplot(1, 2, 1);
imagesc(Im);
axis image;
colormap gray;

subplot(1, 2, 2);
level = il_otsu_threshold(Im, 256);
tic;
BW = il_binary(Im, level);
toc;
imagesc(BW);
axis image;
colormap gray;
