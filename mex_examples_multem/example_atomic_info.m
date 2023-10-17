% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>
clear;clc;
addpath(['..', filesep, 'matlab_functions'])
addpath(['..', filesep, 'crystalline_materials'])
addpath(['..', filesep, 'mex_bin'])

% S = ilc_atomic_info(pot_parm_typ) get atomic properties
% 
% Pot_Par is an integer number which is link to the atomic potential parameterization

% 1: Doyle_0_4
% 2: Peng_0_4
% 3: Peng_0_12
% 4: Kirkland_0_12
% 5: Weickenmeier_0_12
% 6: Lobato_0_12
%
% outputs:
% name              % atom name
% Z                 % atomic number
% A                 % mass number
% m                 % atomic mass
% rn                % experimental nuclear radius (Angs.)
% ra                % experimental atomic radius (Angs.)
% eels_maj_edg      % major eels edges
% eels_min_edg		% minor eels edges
% feg               % electron scattering factor coefficients
% fxg               % x-ray scattering factor coefficients
% pr                % electron density coefficients
% vr                % potential coefficients
% vzp               % projected potential coefficients
% 
% Copyright 2021 Ivan Lobato <Ivanlh20@gmail.com>

clear;clc;
pot_parm_typ = 6;
tic;
S = ilc_atomic_info(pot_parm_typ);
toc;

k = 79;cc = '-+r';
S(k).coef(1).feg
S(k).coef(1).fxg
S(k).coef(1).pr
S(k).coef(1).vr
S(k).coef(1).vzp

S

figure(1); clf;
Z = (1:103)';
y = zeros(size(Z));
for i = 1:5
    hold on;
    subplot(2, 3, i);
    for j = 1:103
        y(j) = S(j).coef(1).feg.cnl(i);
    end
    plot(Z, y, '-*r');
    set(gca, 'FontSize', 12, 'LineWidth', 1, 'PlotBoxAspectRatio', [1.25 1 1]);
    title('Non-Lineal coeeficients');
    ylabel(strcat('cnl[', num2str(i), ']'), 'FontSize', 14);
    xlabel('Z', 'FontSize', 12);
    xlim([1 103]);
    legend(strcat('cnl[', num2str(i), ']'));
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
end