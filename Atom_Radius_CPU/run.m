clear all;
clc;
Dim = 3; Vrl = 0.015;
z = 1:103;
tic;
r = getAtomRadius(6, Dim, Vrl);
toc;
figure(1); clf;
plot(z, sqrt(3)*r(:, 1), '-*r', z, r(:, 2), '-*b', z, r(:, 3), '-*k');
% plot(z, r(:, 1), '-*r', z, r(:, 2), '-*b');
% hold on;
% r = AtomicSize(5, Dim, Vrl);
% plot(z, 3.0*r(:, 1), '-*r', z, r(:, 2), '-*b');

a =0.7; b = 4;
c0 = 0.4378; c1 = 1/c0; c2 = 1-c0; c3 = 1/c2;
[exp(a+b), 1+a+b, c0*exp(c1*a)+c2*exp(c3*b)]