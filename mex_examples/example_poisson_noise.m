clc; clear all;

A = double(rgb2gray(imread('ngc6543a.jpg')));

tic;
B = poissrnd(A);
toc;

SNR = 8;
tic;
[An, k] = get_poisson_noise(A, SNR);
toc;
std(A(:))
k
figure(1);
subplot(1, 2, 1);
imagesc(A);
axis image;
axis image;
subplot(1, 2, 2);
imagesc(B);
axis image;
colormap gray;
axis image;