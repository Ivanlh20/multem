clear all;
clc;
na = 1; nb = 1; nc = 10; pp = 6; ncu = 4; sigma = 0.084; gpu = 0;
Dim = 110; Seed = 1983; iConfFP = 0; nx = 2048; ny = 2048;
[Atomsi, lx, ly, lz, a, b, c, dz] = Au001Crystal(na, nb, nc, ncu, sigma);
Atomsi = [5 5 0 79 0.084 1.0];
 lx = 10; ly = 10; dz = 0.25;
[Atoms, Slice] = getSliceSpecimen(Atomsi, lx, ly, dz, iConfFP, Dim, Seed);
[nSlice, ~] = size(Slice);
for iSlice = 1:nSlice
    [V0, V1] = getProjPotential(Atomsi, gpu, nx, ny, lx, ly, dz, iConfFP, Dim, Seed, iSlice); 
    QuadType = 0; nQuad = 128; PotPar = 6; cPotf = 47.877645145863056;
    [xi, wi] = getQuadrature(QuadType, nQuad);
    Z = 79; Dim = 3; RMSPro = 0.0; sigma = sqrt(RMSPro^2/3); IntTyp = 0;
    S = getAtomTypes(PotPar);

    Rmin = lx/nx; Rmax = S(Z).Rmax; nR = 128;
    dlnR = log(Rmax/Rmin)/(nR-1); Rl = Rmin*exp((0:1:(nR-1))*dlnR);
    V0l = zeros(size(Rl));
    z0 = Slice(iSlice, 1); ze = Slice(iSlice, 2);
    if((z0<0)&&(0<ze))
        a = -0.5*z0;
        b = 0.5*z0;
    else
        a = 0.5*(ze-z0); 
        b = 0.5*(ze+z0);    
    end;
    as = 0.5*ze; 
    bs = 0.5*ze;
    zi = a*xi+b;
    zis = as*xi+bs;
    for i=1:nR
        ri = sqrt(zi.^2+Rl(i)^2);
        [Vi, ~] = getPotential(PotPar, Z, sigma, IntTyp, Dim, ri);
        V0l(i) = a*sum(Vi.*wi)/cPotf;
        V1l(i) = a*sum((zi-z0).*Vi.*wi)/((ze-z0)*cPotf);
        if((z0<0)&&(0<ze))
            ri = sqrt(zis.^2+Rl(i)^2);
            [Vi, ~] = getPotential(PotPar, Z, sigma, IntTyp, Dim, ri);
            V0l(i) = V0l(i) + as*sum(Vi.*wi)/cPotf;
            V1l(i) = V1l(i) + as*sum((zis-z0).*Vi.*wi)/((ze-z0)*cPotf);
        end;
    end;

    dRx = lx/nx; dRy = ly/ny;
    [Rx, Ry] = meshgrid((0:1:(nx-1))*dRx-Atomsi(1,1), (0:1:(ny-1))*dRy-Atomsi(1,2));
    R = sqrt(Rx.^2+Ry.^2); R(R<Rmin) = Rmin;
    V0n = zeros(size(R)); V1n = zeros(size(R));
    ii = find(R<Rmax);
    V0n(ii) = spline(Rl,V0l, R(ii))-V0l(end);
    V1n(ii) = spline(Rl,V1l, R(ii))-V1l(end);

    num2str([sum(abs(V0(:)-V0n(:))), sum(abs(V1(:)-V1n(:)))], 5)
    
    figure(1); clf;
    subplot(1, 2 ,1);
    plot(V0(ny/2+1,:), '-r');
    hold on;
    plot(V0n(ny/2+1,:), '-b');
    subplot(1, 2 ,2);
    plot(V1(ny/2+1,:), '-r');
    hold on;
    plot(V1n(ny/2+1,:), '-b');
    pause(0.25);
end;



[z0, ze]