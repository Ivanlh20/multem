clc; clear all;
L = 20; n = 2048; nh = n/2; dr = L/n; rl = (-nh:1:(nh-1))*dr; typ = 0;
sigma = 2.0; x0 = 0; y0 = 0; rlh = (0:1:(nh-1))*dr;
alpha = 0.5/sigma^2;
[Rx, Ry] = meshgrid(rl, rl); R = sqrt(Rx.^2+Ry.^2);
G2d = exp(-alpha*R.^2);
ya = exp(-alpha*rlh.^2);

rlhn = (0:1.0:(nh-1))*dr;
tic;
[x, f, cf] = il_radial_distribution_2D(R, G2d, rlhn, typ);
toc;
figure(1);
subplot(1, 2,1 );
imagesc(rl, rl, G2d);
colormap gray;
axis image;
subplot(1, 2, 2);
plot(rlh, ya, '-r', x, f, '-b');