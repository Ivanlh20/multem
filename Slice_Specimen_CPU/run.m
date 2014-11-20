clear all;
clc;
na = 5; nb = 5; nc = 100; pp = 6; ncu = 8; sigma = 0;
tic;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);
toc;
Dim = 111; Seed = 1983; inFP = 0;
tic;
[Atoms, Slice] = getSliceSpecimen(Atomsi, lx, ly, dz, inFP, Dim, Seed);
toc;

length(Atoms)

[nSlice, ~] = size(Slice);

for i = 1:nSlice
    figure(1); clf;
    i1 = Slice(i, 7); i2 = Slice(i, 8); ii = i1:1:i2;
    plot3(Atoms(:, 1), Atoms(:, 2), Atoms(:, 3), '+k');
    hold on;
    plot3(Atoms(ii, 1), Atoms(ii, 2), Atoms(ii, 3), '*r');
    axis equal;
    view([1 0 0]);
end;
[length(Atomsi), length(Atoms)]
[lx, ly, lz]