% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
addpath([fileparts(pwd) filesep 'mex_bin'])
addpath([fileparts(pwd) filesep 'crystalline_materials'])
addpath([fileparts(pwd) filesep 'matlab_functions'])

Z = 79;
occ = 1;
region = 0;
charge = 0;

rmin = 1e-02; rmax = 5.0; nr = 512;
dlnr = log(rmax/rmin)/(nr-1); r = rmin*exp((0:1:(nr-1))*dlnr);

tic;
[f1, df1] = ilc_vr(1, Z, charge, r);
[f2, df2] = ilc_vr(2, Z, charge, r);
[f3, df3] = ilc_vr(3, Z, charge, r);
[f4, df4] = ilc_vr(4, Z, charge, r);
[f5, df5] = ilc_vr(5, Z, charge, r);
[f6, df6] = ilc_vr(6, Z, charge, r);
toc;

figure(1); clf;

subplot(1, 2, 1);
hold on;
plot(r, f1, '-k', r, f2, '-b', r, f3, '-c', r, f4, '-m', r, f5, '-r', r, f6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Atomic potential');
ylabel('$\displaystyle V(r)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12);
xlim([0 rmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

subplot(1, 2, 2);
hold on;
plot(r, df1, '-k', r, df2, '-b', r, df3, '-c', r, df4, '-m', r, df5, '-r', r, df6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Derivative of the Atomic potential');
ylabel('$\displaystyle \frac{d V(r)}{dr}$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12);
xlim([0 rmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);