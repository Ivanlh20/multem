close all; clear; clc;
%%
tol         =   1e-12;
a           =   1;
%% func1
rho         =   1;
v           =   0;
tic;
I           =   HankelTransform(@(x)func1(x),rho,v,a,tol);
toc;
% correct answer    :   -1.000000000000000
fprintf('I1\t=\t%0.15f\n',I);
%% func2
rho         =   1;
v           =   1;
tic;
I           =   HankelTransform(@(x)func2(x),rho,v,a,tol);
toc;
% correct answer    :   +0.421024438240708
fprintf('I2\t=\t%0.15f\n',I);
%% func3
rho         =   1;
v           =   0;
tic;
I           =   HankelTransform(@(x)func3(x),rho,v,a,tol);
toc;
% correct answer    :   +1.000000000000000
fprintf('I3\t=\t%0.15f\n',I);
%% func4
rho         =   1;
v           =   10;
tic;
I           =   HankelTransform(@(x)func4(x),rho,v,a,tol);
toc;
% correct answer    :   +0.098970545308402
fprintf('I4\t=\t%0.15f\n',I);
%% func5
rho         =   200;
v           =   0;
tic;
I           =   HankelTransform(@(x)func5(x),rho,v,a,tol);
toc;
% correct answer    :   -1/rho^3
fprintf('I5\t=\t%0.15f\n',I);
%% Test Functions
function[y]=func1(x)
y               =   (x.^2);
end
function[y]=func2(x)
y               =   0.5*log(1+x.^2);
end
function[y]=func3(x)
y               =   (1-exp(-x))./(x.*log(1+sqrt(2)));
end
function[y]=func4(x)
y               =   x./(1+x.^2);
end
function[y]=func5(x)
y               =   x.^2;
end
%%