clear all;
clc;
Z = 79;

gmin = 0; gmax = 12; ng = 512;
dg = (gmax-gmin)/(ng-1); 
g = gmin:dg:gmax;

[f1, df1] = getfeg(1, Z, g);
[f2, df2] = getfeg(2, Z, g);
[f3, df3] = getfeg(3, Z, g);
[f4, df4] = getfeg(4, Z, g);
[f5, df5] = getfeg(5, Z, g);
[f6, df6] = getfeg(6, Z, g);

figure(1);
subplot(1, 2, 1);
plot(g, f3, '-b', g, f4, '-c', g, f5, '-k', g, f6, '-r');
% plot(g, f6/f6(1), '-r');
% plot(g, f1, '-k', g, f2, '-y', g, f3, '-c', g, f4, '-b', g, f5, '-m', g, f6, '-r');
subplot(1, 2, 2);
plot(g, df1, '-k', g, df2, '-y', g, df3, '-c', g, df4, '-b', g, df5, '-m', g, df6, '-r');
% plot(g, f3, '-k', g, f4, '-b', g, f5, '-c', g, f6, '-r');
% xlim([0 200])
% ylim([min(f6) max(f6)]);