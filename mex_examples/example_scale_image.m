clear all; clc;

Im = double(rgb2gray(imread('peppers.png')));
shrink = 0.5;

figure(1) 
subplot(1, 2, 1);
imagesc(imresize(Im, shrink));
axis image;
colormap gray;

subplot(1, 2, 2);
tic;
Im_d = il_scale_image(Im, shrink);
toc;
imagesc(Im_d);
axis image;
colormap gray;
