clear all; clc;
na = 3; nb = 3; nc = 5; pp = 6; ncu = 8; sigma = 0;
tic;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);
toc;
Dim = 111; Seed = 1983; iConfFP = 0;

tic;
% get specimen slicing
[Atoms, Slice] = get_SliceSpecimen_CPU(Atomsi, lx, ly, dz, iConfFP, Dim, Seed);
toc;
[nAtoms,~] = size(Atoms); [nSlice, ~] = size(Slice);

for i = 1:nSlice
    figure(1); clf;
    i1 = Slice(i, 7); i2 = Slice(i, 8); ii = i1:1:i2;
    plot3(Atoms(:, 1), Atoms(:, 2), Atoms(:, 3), '*k', Atoms(ii, 1), Atoms(ii, 2), Atoms(ii, 3), '*r');
    set(gca,'FontSize',12,'LineWidth',1,'PlotBoxAspectRatio',[1.25 1 1]);
    title('Atomic positions');
    ylabel('y','FontSize',14);
    xlabel('x','FontSize',12);
    axis equal;
    view([1 0 0]);
    pause(0.1);
end;

[length(Atomsi), length(Atoms)]
[lx, ly, lz]