clc; clear all;
L = 20; n = 2048; nh = n/2; dr = L/n; rl = (-nh:1:(nh-1))*dr; typ = 0;
sigma = 2.0; x0 = 0; y0 = 0; rlh = (0:1:(nh-1))*dr;
alpha = 0.5/sigma^2;
[Rx, Ry] = meshgrid(rl, rl); R = sqrt(Rx.^2+Ry.^2);
G2d = exp(-alpha*R.^2);
fG2d = abs(fft2(ifftshift(G2d)));
radius = get_FFT_InformationLimit_2D_CPU(fG2d, 0)