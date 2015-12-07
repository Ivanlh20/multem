clear all; clc;

Im = double(rgb2gray(imread('peppers.png')));
level = 0.5*(max(Im(:))+min(Im(:)));
BW = Im;
BW(BW<level) = 0;
BW(BW>=level) = 1;

figure(1) 
subplot(1, 3, 1);
imagesc(Im);
axis image;
colormap gray;

subplot(1, 3, 2);
imagesc(BW);
axis image;
colormap gray;

level = il_otsu_threshold(Im, 256);
BW = Im;
BW(BW<level) = 0;
BW(BW>=level) = 1;
subplot(1, 3, 3);
imagesc(BW);
axis image;
colormap gray;
% nbins = 100;
% 
% figure(1); clf;
% subplot(1, 2, 1)
% tic;
% hist(Im,nbins);
% toc;
% 
% tic;
% [x, y] = il_histogram(Im, nbins);
% toc;
% hold on;
% plot(x, y, '-*r');

