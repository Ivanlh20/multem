clear all; clc;
Z = 6; dg = 2;
% Load Kirkland X-ray and electron scattering tabulated data
tic;
[g, fxg, feg] = get_fxegTabData_CPU(Z, dg);
toc;

figure(1); clf;

subplot(1, 2, 1);
plot(g, fxg, '-+r');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Tabulated X-ray Scattering factor');
ylabel('$\displaystyle f_x(g)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{g}$','interpreter','latex','FontSize',12,'FontWeight','bold');
xlim([0 g(end)]);
legend('Kirkland [0-12]');

subplot(1, 2, 2);
plot(g, feg, '-+r');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Tabulated electron Scattering factor');
ylabel('$\displaystyle f_e(g)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{g}$','interpreter','latex','FontSize',12);
xlim([0 g(end)]);
legend('Kirkland [0-12]');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);