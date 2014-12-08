clear all; clc;

na = 1; nb = 1; nc = 5; pp = 6; ncu = 2; sigma = 0.084;
InProjPotential.gpu = 0;                % Gpu card
InProjPotential.iConfFP = 0;            % Frozen phonon configuration
InProjPotential.DimFP = 110;            % Dimensions phonon configuration
InProjPotential.SeedFP = 1983;          % Frozen phonon random seed
InProjPotential.PotPar = 6;             % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
InProjPotential.nx = 2048; InProjPotential.ny = 2048;
[InProjPotential.Atoms, InProjPotential.lx, InProjPotential.ly, lz, a, b, c, InProjPotential.dz] = Au001Crystal(na, nb, nc, ncu, sigma);
InProjPotential.Atoms = [4.0 4.0 0 79 sigma 1.0];
InProjPotential.lx = 8.0; InProjPotential.ly = 8.0; InProjPotential.dz = 0.5;
[Atoms, Slice] = getSliceSpecimen(InProjPotential.Atoms, InProjPotential.lx, InProjPotential.ly, InProjPotential.dz, InProjPotential.iConfFP, InProjPotential.DimFP, InProjPotential.SeedFP);
[nAtoms,~] = size(Atoms); [nSlice, ~] = size(Slice);
for iSlice = 1:nSlice
    InProjPotential.iSlice = iSlice;
    tic;
    clear getProjPotential;
    [V0, V1] = getProjPotential(InProjPotential);
    toc;
    figure(1);
    subplot(1, 2, 1);    
    imagesc(V0);
    colormap gray;
    axis image;
    subplot(1, 2, 2);    
    imagesc(V1);
    colormap gray;
    axis image;    
    pause(0.10);
end;
