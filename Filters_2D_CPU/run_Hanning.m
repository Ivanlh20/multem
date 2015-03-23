clc; clear all;
nx = 1024; ny = 1024;
dx = 1; dy = 1; k = 0.1; shift = 1;
f = get_Hanning_Filter_2D_CPU(ny, nx, dy, dx, k, shift);
figure(1);
imagesc(f);
colormap gray;
axis image;