clear all; clc;
Dim = 3; Vrl = 0.015; PotPar = 6;
z = 1:103;
tic;
r = getAtomRadius(PotPar, Dim, Vrl);
toc;

figure(1); clf;

plot(z, sqrt(3)*r(:, 1), '-*r', z, r(:, 2), '-*b', z, r(:, 3), '-*k');
set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
title('Atomic radius');
ylabel('$\mathrm{radius}$','interpreter','latex','FontSize',14);
xlabel('$\mathbf{r}$','interpreter','latex','FontSize',12);
xlim([1 103]);
legend('rms', 'Cut-off', 'Experimental');

set(gcf,'units','normalized','outerposition',[0 0 1 1]);