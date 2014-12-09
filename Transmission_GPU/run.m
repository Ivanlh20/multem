clear all; clc;

InTransmission.gpu = 0;             % Gpu card
InTransmission.iConfFP = 0;         % Frozen phonon configuration
InTransmission.DimFP = 110;         % Dimensions phonon configuration
InTransmission.SeedFP = 1983;       % Frozen phonon random seed
InTransmission.PotPar = 6;          % Parameterization of the potential 1: Doyle(0-4), 2: Peng(0-4), 3: peng(0-12), 4: Kirkland(0-12), 5:Weickenmeier(0-12) adn 6: Lobato(0-12)
InTransmission.ApproxModel = 1;     % 1: MS, 2: PA, 3:POA, 4:WPOA
InTransmission.BWL = 1;             % 1: true, 2: false
InTransmission.FastCal = 1;         % 1: normal mode(low memory consumption), 2: fast calculation(high memory consumption)

InTransmission.E0 = 300;
InTransmission.theta = 0.0; InTransmission.phi = 0.0; % Till ilumination (degrees)
InTransmission.nx = 2048; InTransmission.ny = 2048;

na = 3; nb = 3; nc = 5; pp = 6; ncu = 4; sigma = 0.084;
[InTransmission.Atoms, InTransmission.lx, InTransmission.ly, lz, a, b, c, InTransmission.dz] = Au001Crystal(na, nb, nc, ncu, sigma);
% Atomsi = [ 2.0 2.0 0 79 0.084 1.0];
%  lx = 4; ly = 4.0; dz = 0.5;
[Atoms, Slice] = getSliceSpecimen(InTransmission.Atoms, InTransmission.lx, InTransmission.ly, InTransmission.dz, InTransmission.iConfFP, InTransmission.DimFP, InTransmission.SeedFP);
[nAtoms,~] = size(Atoms); [nSlice, ~] = size(Slice);
for iSlice = 1:nSlice
    InTransmission.iSlice = iSlice;
    tic;
    clear getTransmission;
    Trans = getTransmission(InTransmission);
    toc;
    figure(1);
    subplot(1, 3, 1);    
    imagesc(real(Trans));
    colormap gray;
    axis image;
    subplot(1, 3, 2);    
    imagesc(imag(Trans));
    colormap gray;
    axis image;   
    subplot(1, 3, 3);    
    imagesc(abs(Trans));
    colormap gray;
    axis image;   
    num2str([iSlice, min(abs(Trans(:))), max(abs(Trans(:))), sum(abs(Trans(:)))/(InTransmission.nx*InTransmission.ny)], 10)
    pause(0.10);
end;