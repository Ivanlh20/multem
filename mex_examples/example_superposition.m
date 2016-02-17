clc; clear all;

%read peaks
nx = 1024;
ny = 1024;
dx = 1;
dy = 1;

%% read points and scale it
xy = importdata('peaks.mat');
xy(:, 1) = xy(:, 1)*nx*dx/2048;
xy(:, 2) = xy(:, 2)*ny*dy/2048;

n = size(xy, 1);
height = 1;
sigma = 2;
ff_sigma = 4;
data = [xy, height*ones(n, 1), sigma*ones(n, 1)]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Rx, Ry] = meshgrid((0:1:(nx-1))*dx, (0:1:(ny-1))*dy);
Im_1 = zeros(ny, nx);
R2_max = (ff_sigma*sigma)^2;
tic;
for i = 1:n;
    R2 = (Rx-data(i, 1)).^2 + (Ry-data(i, 2)).^2;
    ii = find(R2<R2_max);
    Im_1(ii) = Im_1(ii) + data(i, 3)*exp(-0.5*R2(ii)/data(i, 4)^2);
end;
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C++%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
superposition.precision = 2;                     % eP_Float = 1, eP_double = 2
superposition.device = 1;                        % eD_CPU = 1, eD_GPU = 2
superposition.cpu_nthread = 4; 
superposition.gpu_device = 0;
superposition.gpu_nstream = 8;

superposition.ff_sigma = ff_sigma;
superposition.nx = nx;
superposition.ny = ny;
superposition.dx = 1;
superposition.dy = 1;
superposition.data =  data;

tic;
Im_2 = il_superposition(superposition);
toc;

figure(1);
subplot(1, 2, 1);
imagesc(Im_1);
axis image;
subplot(1, 2, 2);
imagesc(Im_2);
axis image;

figure(2);
imagesc(Im_1-Im_2);