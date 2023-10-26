% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

Z = 49;
occ = 1;
region = 0;
charge = 0;

gmin = 0; gmax = 20; ng = 512;
dg = (gmax-gmin)/(ng-1); 
g = gmin:dg:gmax;

[f1, df1] = ilc_fxg(1, Z, charge, g);
[f2, df2] = ilc_fxg(2, Z, charge, g);
[f3, df3] = ilc_fxg(3, Z, charge, g);
[f4, df4] = ilc_fxg(4, Z, charge, g);
[f5, df5] = ilc_fxg(5, Z, charge, g);
[f6, df6] = ilc_fxg(6, Z, charge, g);

figure(1); clf;

subplot(1, 2, 1);
plot(g, f1, '-k', g, f2, '-b', g, f3, '-c', g, f4, '-m', g, f5, '-r', g, f6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('X-ray Scattering factor');
ylabel('$\displaystyle f_x(g)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{g}$','interpreter','latex','FontSize',12);
xlim([0 gmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

subplot(1, 2, 2);
plot(g, df1, '-k', g, df2, '-b', g, df3, '-c', g, df4, '-m', g, df5, '-r', g, df6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Derivative of the X-ray Scattering factor');
ylabel('$\displaystyle \frac{d f_x(g)}{dg}$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{g}$','interpreter','latex','FontSize',12);
xlim([0 gmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);