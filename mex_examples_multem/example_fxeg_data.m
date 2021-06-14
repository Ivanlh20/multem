% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath([ fileparts(pwd), filesep, 'mex_bin'])

Z = 24; dg = 1;
% Load Kirkland X-ray and electron scattering tabulated data
tic;
[g, fxg, feg] = ilc_fxeg_data(Z, dg);
toc;

a = [3.028317848436913e+00	-9.553939330814323e+01	9.617615623521981e+01	3.597773159579877e-02	3.914928907184223e-04];
b = [8.359115043146318e+00	1.802637902641032e+00	1.775094889570472e+00	5.481444121431137e-02	3.998289689160348e-03];
a_0 = 0.52917721077817892;

feg_p = sum(a.*(2+b.*g.^2)./(1+b.*g.^2).^2, 2);
% feg = (Z-fxg)./(2*pi^2*a_0*g.^2);

figure(1); clf;
plot(g, feg, '-+r', g, feg_p, '-+b');
set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
title('Tabulated electron Scattering factor');
ylabel('$\displaystyle f_e(g)$', 'interpreter', 'latex', 'FontSize', 14);
xlabel('$\mathbf{g}$', 'interpreter', 'latex', 'FontSize', 12);
xlim([0 g(end)]);
legend('Kirkland [0-12]');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);