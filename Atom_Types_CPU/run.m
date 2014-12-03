%Ivan lobato - 26/06/2012-6.30pm
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PotPar = 6;
tic;
S = getAtomTypes(PotPar);
toc;
k = 6; cc = '-+r';
S(k)

figure(1); clf;
Z = (1:103)';
y = zeros(size(Z));
for i = 1:5
    hold on;
    subplot(2, 3, i);
    for j = 1:103;
        y(j) = S(j).cVr.cnl(i);
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