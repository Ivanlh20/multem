clear all;
clc;
na = 1; nb = 1; nc = 10; pp = 6; ncu = 4; sigma = 0.084; gpu = 0;
Dim = 110; Seed = 1983; iConfFP = 0; nx = 2048; ny = 2048;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);
Atomsi = [ 5 5 0 79 0.084 1.0];
 lx = 10; ly = 10; dz = 0.25;
[Atoms, Slice] = getSliceSpecimen(Atomsi, lx, ly, dz, iConfFP, Dim, Seed);

[nAtoms,~] = size(Atoms);
[nSlice, ~] = size(Slice);
for iSlice = 1:nSlice
%     i1 = Slice(iSlice, 3); i2 = Slice(iSlice, 4); ii1 = i1:1:i2;
%     i1 = Slice(iSlice, 7); i2 = Slice(iSlice, 8); ii2 = i1:1:i2;
%     subplot(1, 4, 2);
%     plot3(Atoms(:, 1), Atoms(:, 2), Atoms(:, 3), '*k', Atoms(ii2, 1), Atoms(ii2, 2), Atoms(ii2, 3), '*r', Atoms(ii1, 1), Atoms(ii1, 2), Atoms(ii1, 3), '*b');
%     view ([1 0 0]);
%     axis equal;
    
    tic;
    [V0, V1] = getProjPotential(Atomsi, gpu, nx, ny, lx, ly, dz, iConfFP, Dim, Seed, iSlice);
    toc;
    figure(1);
    subplot(1, 2, 1);    
    imagesc(V0);
    colormap hot;
    axis image;
    subplot(1, 2, 2);    
    imagesc(V1);
    colormap hot;
    axis image;    
    pause(0.25);
end;
