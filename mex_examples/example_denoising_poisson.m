clc; clear all;
RGB = imread('saturn.png');
I = 0.1*double(rgb2gray(RGB));
[ny, nx] = size(I);

J = poissrnd(I);
% J = I + 25*randn(ny, nx);

[min(I(:)), max(I(:))]

figure(1);
subplot(2, 2, 1);
imagesc(I);
axis image;
colormap gray;

subplot(2, 2, 2);
imagesc(J);
axis image;
colormap gray;

tic;
K = wiener2(J,[7 7]);
toc;
[min(K(:)), max(K(:))]
subplot(2, 2, 3);
imagesc(K);
axis image;
colormap gray;

tic;
L = il_denoising_poisson(J, 3, 2);
% L = il_anscombe_forward(J);
% L = il_filter_wiener(L, 3);
% L = il_filter_median(L, 2);
% L = il_anscombe_inverse(L);
toc;
[min(L(:)), max(L(:))]

subplot(2, 2, 4);
imagesc(L);
axis image;
colormap gray;

figure(2);
x = 1:nx;
plot(x, J(700, :), '-k', x, K(700, :), '-b', x, L(700, :), '-c', x, I(700, :), '-r');
