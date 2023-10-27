% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;

E_0 = [60, 80, 100, 120, 200, 300];
sigma = ilc_elec_interact_parm_kva(E_0); % radians/(kV - Angs)
tf_exp_factor = ilc_transf_exp_factor(E_0); % radians/(V - Angs)

f = 1e-3

figure(1); clf;
plot(E_0, tf_exp_factor, '-*r', E_0, sigma*f, '-*b');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
% title('Atomic radius');
% ylabel('$\mathrm{radius}$', 'interpreter', 'latex', 'FontSize', 14);
% xlabel('$\mathbf{r}$', 'interpreter', 'latex', 'FontSize', 12);
% axis([1 103 0 1.1*max(r(:))]);
% legend('rms', 'Cut-off', 'Experimental');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);