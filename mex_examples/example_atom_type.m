% S = il_atom_type(PotPar) get atomic properties
% 
% Pot_Par is an integer number which is link to the atomic potential parameterization
% 
%     1: Doyle_0_4
%     2: Peng_0_4
%     3: Peng_0_12
%     4: Kirkland_0_12
%     5: Weickenmeier_0_12
%     6: Lobato_0_12
% 
% Dim is the dimension with possible values 2 and 3
% 
% Vrl is the atomic potential threshold in eV
% 
% r(:, 1) = integral {r^2 V(r)dr} / integral{V(r)dr}
% r(:, 2) are calculated by solving the equation V(r) = Vrl
% r(:, 3) are the experimental values
% 
% Copyright 2016 Ivan Lobato <Ivanlh20@gmail.com>

clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PotPar = 6;
tic;
S = il_atom_type(PotPar);
toc;
k = 79; cc = '-+r';
S(k).coef(1).feg
S(k).coef(1).fxg
S(k).coef(1).Pr
S(k).coef(1).Vr
S(k).coef(1).VR

figure(1); clf;
Z = (1:103)';
y = zeros(size(Z));
for i = 1:5
    hold on;
    subplot(2, 3, i);
    for j = 1:103;
        y(j) = S(j).coef(1).feg.cnl(i);
    end;
    plot(Z, y, '-*r');
    set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
    title('Non-Lineal coeeficients');
    ylabel(strcat('cnl[',num2str(i),']'), 'FontSize',14);
    xlabel('Z','FontSize',12);
    xlim([1 103]);
    legend(strcat('cnl[',num2str(i),']'));  
    set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
end;