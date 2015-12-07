clc; clear all;

% A = double(rgb2gray(imread('ngc6543a.jpg')));
nx = 1024;
lx = 20;
lxh = lx/2;
dx = lx/nx;
d = 4;

[Rx, Ry] = meshgrid(0:dx:lx);
x =[lxh-d, lxh; lxh+d, lxh];
x =[lxh, lxh];
sigma = 2;
A0 = zeros(size(Rx));
for i = 1:1
    A0 = A0 + exp(-0.5*((Rx-x(i, 1)).^2+(Ry-x(i, 2)).^2)/sigma^2);
end;

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

% dRx = 1;
% dRy = 1;
% sigma = 5;
% tic;
% Ab = il_convolution_gaussian(A, dRx, dRy, sigma);
% toc;
% 
% figure(1);
% subplot(1, 2, 1);
% imagesc(A);
% axis image;
% axis image;
% subplot(1, 2, 2);
% imagesc(Ab);
% axis image;
% colormap gray;
% axis image;