clear all; clc;

Z = 79;

Rmin = 1e-02; Rmax = 5.0; nR = 512;
dlnR = log(Rmax/Rmin)/(nR-1); R = Rmin*exp((0:1:(nR-1))*dlnR);

tic;
[f1, df1] = get_Vp(1, Z, R);
[f2, df2] = get_Vp(2, Z, R);
[f3, df3] = get_Vp(3, Z, R);
[f4, df4] = get_Vp(4, Z, R);
[f5, df5] = get_Vp(5, Z, R);
[f6, df6] = get_Vp(6, Z, R);
toc;

figure(1); clf;

subplot(1, 2, 1);
hold on;
plot(R, f1, '-k', R, f2, '-b', R, f3, '-c', R, f4, '-m', R, f5, '-r', R, f6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Projected Atomic potential');
ylabel('$\displaystyle V(R)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{R}$','interpreter','latex','FontSize',12);
xlim([0 Rmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

subplot(1, 2, 2);
hold on;
plot(R, df1, '-k', R, df2, '-b', R, df3, '-c', R, df4, '-m', R, df5, '-r', R, df6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Derivative of the Atomic potential');
ylabel('$\displaystyle \frac{d V(R)}{dr}$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{R}$','interpreter','latex','FontSize',12);
xlim([0 Rmax]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);