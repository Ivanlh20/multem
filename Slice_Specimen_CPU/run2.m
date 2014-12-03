clear all;
clc;
na = 10; nb = 10; nc = 5; PotPar = 6; ncu = 4; sigma = 0.085;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);

Dim = 111; Seed = 1983; iConfFP = 2;
tic;
% get specimen slicing
[Atoms, Slice] = getSliceSpecimen(Atomsi, lx, ly, dz, iConfFP, Dim, Seed);
toc;
[nAtoms,~] = size(Atoms); [nSlice, ~] = size(Slice);
S = getAtomTypes(PotPar);
z0 = min(Atoms(:, 3))-S(Atoms(1,4)).Rmax; 
ze = max(Atoms(:, 3))+S(Atoms(end,4)).Rmax;

xy = zeros(nAtoms, 2);
xy(:, 2) = Atoms(:, 3);
figure(1); clf;
plot(xy(:, 1), xy(:, 2), '*b');
hold on;
plot([-1 1], [z0 z0], '-k','LineWidth',2);
hold on;
plot([-1 1], [ze ze], '-k','LineWidth',2);
for i = 1:nSlice
    hold on;
    plot([-1 1], [Slice(i, 1) Slice(i, 1)], '-r');
    set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[0.75 1 1]);
    title('Atomic positions');
    ylabel('z','FontSize',14);
    xlabel('x','FontSize',12);
end;
hold on;
plot([-1 1], [Slice(i, 2) Slice(i, 2)], '-r');

