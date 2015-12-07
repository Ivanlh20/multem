clear all; clc;

Im = randn(2048^2,1);
nbins = 100;

figure(1); clf;
subplot(1, 2, 1)
tic;
hist(Im,nbins);
toc;

tic;
[x, y] = il_histogram(Im, nbins);
toc;
hold on;
plot(x, y, '-.r');

