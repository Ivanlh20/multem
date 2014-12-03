clear all; clc;
a0 = 0.52917721077817892;
Z = 79;
rmin = 1e-05; rmax = 0.1; nr = 512;
dlnr = log(rmax/rmin)/(nr-1); r = rmin*exp((0:1:(nr-1))*dlnr);

[f1, df1] = getrho(1, Z, r); 
[f2, df2] = getrho(2, Z, r);
[f3, df3] = getrho(3, Z, r);
[f4, df4] = getrho(4, Z, r);
[f5, df5] = getrho(5, Z, r);
[f6, df6] = getrho(6, Z, r);

figure(1);
subplot(1, 2, 1);
plot(r, f1, '-k', r, f2, '-b', r, f3, '-c', r, f4, '-m', r, f5, '-r', r, f6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Electron density');
ylabel('$\displaystyle \rho(r)$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12);
xlim([0 rmax]); ylim([0 max(f6)]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

subplot(1, 2, 2);
plot(r, df1, '-k', r, df2, '-b', r, df3, '-c', r, df4, '-m', r, df5, '-r', r, df6, '-g');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Derivative of the Electron density');
ylabel('$\displaystyle \frac{d \rho(r)}{d r}$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12);
xlim([0 rmax]); ylim([min(df6) max(df6)]);
legend('Doyle [0-4]', 'Peng [0-4]', 'Peng [0-12]', 'Kirkland [0-12]', 'Weickenmeier [0-12]', 'Lobato [0-12]');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);