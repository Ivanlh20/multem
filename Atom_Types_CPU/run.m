%Ivan lobato - 26/06/2012-6.30pm
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PotPar = 6;
tic;
S = getAtomTypes(PotPar);
toc;
k = 79; cc = '-+r';
S(k)
S(k).cfeg
S(k).cfxg
S(k).cPr
S(k).cVr
S(k).cVR

figure(1); clf;
Z = (1:103)';
for i = 1:5
    hold on;
    subplot(2, 3, i);
    for j = 1:103;
        y(j) = S(j).cVr.cnl(i);
    end;
    plot(Z, y, '-*r');
    set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.5 1 1]);
end;