clc; clear all;
RGB = imread('saturn.png');
I = double(rgb2gray(RGB));
[ny, nx] = size(I);
J = I + 15*randn(ny, nx);

figure(1);
subplot(1, 4, 1);
imagesc(J);
axis image;
colormap gray;

tic;
L = il_filter_wiener(J,2);
toc;
subplot(1, 4, 2);
imagesc(L);
axis image;
colormap gray;

tic;
Lr = il_filter_wiener_by_row(J, 2);
toc;
subplot(1, 4, 3);
imagesc(Lr);
axis image;
colormap gray;

tic;
L = J;
for iy =1:ny
    L(iy,:) = il_filter_wiener(J(iy,:), 2);
end;
toc;
subplot(1, 4, 4);
imagesc(L);
axis image;
colormap gray;

sum(abs(L(:)-Lr(:)))