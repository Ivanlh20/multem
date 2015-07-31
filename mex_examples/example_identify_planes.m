clear all; clc;
% na = 5; nb = 5; nc = 5; ncu = 4; PotPar = 6; DWAu3D = 0.6373; rmsAu3D = sqrt(DWAu3D/(8*pi^2)); sigma = rmsAu3D;
% [Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);   

na = 8; nb = 8; nc = 4; ncu = 4; sigma = 0.0;
[Atomsi, lx, ly, lz, a, b, c, dz] = SrTiO3001Crystal(na, nb, nc, ncu, sigma);
Atomsi(Atomsi(:, 4)==8, 5) = sqrt(0.009275); Atomsi(Atomsi(:, 4)==22, 5) = sqrt(0.00557); 
Atomsi(Atomsi(:, 4)==38, 5) = sqrt(0.00787);
 
Dim = 111; Seed = 300183; iConfFP = 2; PotPar = 6;
tic;
% get specimen slicing
[Atoms, Slice] = get_specimen_slicing(Atomsi, lx, ly, dz, iConfFP, Dim, Seed);
toc;
[nAtoms,~] = size(Atoms); [nslice, ~] = size(Slice);
S = get_atom_type(PotPar);
z0 = min(Atoms(:, 3))-S(Atoms(1,4)).Rmax; 
ze = max(Atoms(:, 3))+S(Atoms(end,4)).Rmax;

xy = zeros(nAtoms, 2);
xy(:, 1) = Atoms(:, 1);
xy(:, 2) = Atoms(:, 3);

figure(1); clf;
plot(xy(:, 1), xy(:, 2), 'ok','LineWidth',2);
hold on;
ee = 0e-02;
plot([-2 ly+2], [z0-ee z0-ee], '-k','LineWidth',3);
hold on;
plot([-2 ly+2], [ze+ee ze+ee], '-k','LineWidth',3);
for i = 1:nslice
    hold on;
    plot([-2 ly+2], [Slice(i, 1) Slice(i, 1)], '-r','LineWidth',1);
    hold on;
    plot([-2 ly+2], [Slice(i, 2) Slice(i, 2)], '-c','LineWidth',1);
    set(gca,'FontSize',10,'PlotBoxAspectRatio',[0.75 1 1]);
end;
hold on;
plot([-2 ly+2], [Slice(i, 2) Slice(i, 2)], '-r','LineWidth',1);

tic;
Planes = get_Identify_Planes_CPU(Atomsi, lx, ly, iConfFP, Dim, Seed);
toc;
[nPlanes, ~] = size(Planes);
for i = 1:nPlanes
    hold on;
    plot([-2 ly+2], [Planes(i) Planes(i)], '-b');
    set(gca,'FontSize',10,'PlotBoxAspectRatio',[0.75 1 1]);
end;
[nAtoms, nslice, nPlanes]
% 
