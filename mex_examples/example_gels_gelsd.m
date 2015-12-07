clc; clear all;

t = [0 .3 .8 1.1 1.6 2.3]';
y = [.82 .72 .63 .60 .55 .50]';

A = [ones(size(t)) exp(-t)];
% matlab
tic;
x1 = A\y;
toc;
% mex
tic;
x2 = il_gels(A,y);
toc;
tic;
x3 = il_gelsd(A,y);
toc;

ti = (0:0.1:2.5)';
figure(1);

yi1 = [ones(size(ti)) exp(-ti)]*x1;
yi2 = [ones(size(ti)) exp(-ti)]*x2;
yi3 = [ones(size(ti)) exp(-ti)]*x3;
plot(ti, yi1, '-k', ti, yi2, '-b',  ti, yi2, '-r', t, y,'o')
