clear; clc;
breaks = [0, 1, 2];
c = [2, 1, 1; 3 1 4];
pp = mkpp(breaks, c);
x = 0:1e-4:2;
figure(1); clf;
plot(x, ppval(pp, x), '-r');