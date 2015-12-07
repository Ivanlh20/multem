clc; clear all;
RGB = imread('saturn.png');
I = double(rgb2gray(RGB));
[ny, nx] = size(I);
J = I + 15*randn(ny, nx);

figure(1);
subplot(1, 3, 1);
imagesc(J);
axis image;
colormap gray;

tic;
K = wiener2(J,[5 5]);
toc;
subplot(1, 3, 2);
imagesc(K);
axis image;
colormap gray;

tic;
L = il_filter_mwiener(J,2);
toc;
subplot(1, 3, 3);
imagesc(L);
axis image;
colormap gray;
