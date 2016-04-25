clc; clear all;

% A = double(rgb2gray(imread('ngc6543a.jpg')));
nx = 1024;
ny = 1024;
lx = 20;
ly = 20;
dx = lx/nx;
dy = ly/ny;

dx = lx/nx; dy = ly/ny; sigma = 1; shift = 0;

tic;
Im = il_func_gaussian(nx, ny, dx, dy, sigma, shift);


figure(1);
subplot(1, 2, 1)
imagesc(A0);
axis image;
colormap gray;


dRx = dx;
dRy = dx;
sigma = 2.0;
tic;
Ab = il_convolution_gaussian(A0, dRx, dRy, sigma);
toc;
subplot(1, 2, 2)
imagesc(Ab);
axis image;
colormap gray;