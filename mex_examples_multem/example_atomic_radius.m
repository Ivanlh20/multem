% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

% r = ilc_atomic_radius(Pot_Par, dim, Vrl) calculates the atomic radius
% 
% Pot_Par is an integer number which is link to the atomic potential parameterization
% 
%  1: Doyle_0_4
%  2: Peng_0_4
%  3: Peng_0_12
%  4: Kirkland_0_12
%  5: Weickenmeier_0_12
%  6: Lobato_0_12
% 
% dim is the dimension with possible values 2 and 3
% 
% Vrl is the atomic potential threshold in eV
% 
% r(:, 1) = integral {r^2 V(r)dr} / integral{V(r)dr}
% r(:, 2) are calculated by solving the equation V(r) = Vrl
% r(:, 3) are the experimental values
% 
% Copyright 2023 Ivan Lobato <Ivanlh20@gmail.com>

clear; clc;
dim = 3;
Vrl = 0.01;
pot_parm_typ = 6;
z = 1:103;
tic;
r = ilc_atomic_radius(pot_parm_typ, dim, Vrl);
toc;

figure(1); clf;

plot(z, r(:, 1), '-*r', z, r(:, 2), '-*b', z, r(:, 3), '-*k');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
title('Atomic radius');
ylabel('$\mathrm{radius}$', 'interpreter', 'latex', 'FontSize', 14);
xlabel('$\mathbf{r}$', 'interpreter', 'latex', 'FontSize', 12);
axis([1 103 0 1.1*max(r(:))]);
legend('rms', 'Cut-off', 'Experimental');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);