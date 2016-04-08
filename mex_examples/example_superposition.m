clc; clear all;

%read peaks
nx = 1024;
ny = 1024;
dx = 1;
dy = 1;

%% read peaks and scale it
data = importdata('peaks.mat');
data(:, 1) = data(:, 1)*nx*dx/2048;
data(:, 2) = data(:, 2)*ny*dy/2048;
data(:, 4) = data(:, 4)*nx*dx/2048;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C++%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
superposition.precision = 2;                     % eP_Float = 1, eP_double = 2
superposition.device = 1;                        % eD_CPU = 1, eD_GPU = 2
superposition.cpu_nthread = 4; 
superposition.gpu_device = 0;
superposition.gpu_nstream = 8;

superposition.ff_sigma = 4;
superposition.nx = nx;
superposition.ny = ny;
superposition.dx = 1;
superposition.dy = 1;
superposition.data =  data;

tic;
Im = il_superposition(superposition);
toc;

figure(1);
imagesc(Im);
colormap gray;
axis image;