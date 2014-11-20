clear all;
clc;
na = 4; nb = 4; nc = 10; PotPar = 6; ncu = 2; sigma = 0.085;
tic;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);
Atomsi = [5 5 0.0 79 sigma 1.00; 5 5 5.0 79 sigma 1.00];
Atomsi = [5 5 0 79 0.084 1.0];
Atomsi = [5 5 0.0 79 sigma 1.0;5 5 2.5 79 sigma 1.0;5 5 5.0 79 sigma 1.0;5 5 7.5 79 sigma 1.0;5 5 10 79 sigma 1.0];
lx = 10; ly = 10;

na = 20; nb = 14; nc = 40; ncu = 2; DWAu3D = 0.6373; rmsAu3D = sqrt(DWAu3D/(8*pi^2)); sigma = sqrt(rmsAu3D^2/3);
[Atomsi, lx, ly, lz, a, b, c, dz] = Au110Crystal(na, nb, nc, ncu, sigma);   


dz = 4.078/(2*sqrt(2));
toc;
Dim = 001; Seed = 1987; inFP = 0;
tic;
[Atoms, Slice] = getSliceSpecimen(Atomsi, lx, ly, dz, inFP, Dim, Seed);
[nSlice, ~] = size(Slice)
toc;
S = getAtomTypes(PotPar);
z0 = min(Atoms(:, 3))-S(Atoms(1,4)).Rmax;
ze = max(Atoms(:, 3))+S(Atoms(end,4)).Rmax;

[nAtoms,~] = size(Atoms);
[nSlice, ~] = size(Slice);

xy = zeros(nAtoms, 2);
xy(:, 1) = Atoms(:, 1);
xy(:, 2) = Atoms(:, 3);


figure(1); clf;
plot(xy(:, 1), xy(:, 2), '*k');
hold on;
ee = 0e-02;
plot([-2 ly+2], [z0-ee z0-ee], '-k','LineWidth',3);
hold on;
plot([-2 ly+2], [ze+ee ze+ee], '-k','LineWidth',3);
for i = 1:nSlice
    hold on;
    plot([-2 ly+2], [Slice(i, 1) Slice(i, 1)], '-r','LineWidth',1);
    hold on;
    plot([-2 ly+2], [Slice(i, 2) Slice(i, 2)], '-c','LineWidth',1);
    set(gca,'FontSize',10,'PlotBoxAspectRatio',[0.75 1 1]);
end;
hold on;
plot([-2 ly+2], [Slice(i, 2) Slice(i, 2)], '-r','LineWidth',1);
tic;
Planes = getPlanes(Atomsi, lx, ly, inFP, Dim, Seed);
toc;
[nPlanes, ~] = size(Planes);
for i = 1:nPlanes
    hold on;
    plot([-2 ly+2], [Planes(i) Planes(i)], '-b');
    set(gca,'FontSize',10,'PlotBoxAspectRatio',[0.75 1 1]);
end;
[nAtoms, nSlice, nPlanes]
% 
