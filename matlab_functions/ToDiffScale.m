function [DPsi] = ToDiffScale(fPsi)
c = 0.1;
DPsi = log(1+c*abs(fPsi).^2);