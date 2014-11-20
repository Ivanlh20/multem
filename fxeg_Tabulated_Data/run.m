clear all;
clc;
Z = 6;
[g, fxg, feg] = getfxegTabData(Z, 2);
figure(1);
subplot(2, 1, 1);
plot(g, fxg, '+r');

subplot(2, 1, 2);
plot(g, feg, '+r');