function []=ilm_show_slicing(input_multem, inFP, nf)
    [Atoms, Slice] = getSliceSpecimen(input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, input_multem.spec_dz, inFP, input_multem.DimFP, input_multem.SeedFP);
    S = getAtomTypes(input_multem.PotPar);
    z0 = min(Atoms(:, 3))-S(Atoms(1,4)).Rmax;
    ze = max(Atoms(:, 3))+S(Atoms(end,4)).Rmax;
    [nAtoms,~] = size(Atoms);
    [nslice, ~] = size(Slice);

    xy = zeros(nAtoms, 2);
    xy(:, 2) = Atoms(:, 3);
    figure(nf); clf;
    plot(xy(:, 1), xy(:, 2), '*k');
    for i = 1:nslice
        hold on;
        plot([-1 1], [Slice(i, 1) Slice(i, 1)], '-r','LineWidth',1);
        set(gca,'FontSize',10,'LineWidth',2,'PlotBoxAspectRatio',[0.75 1 1]);
    end
    hold on;
    plot([-1 1], [Slice(i, 2) Slice(i, 2)], '-r','LineWidth',1);
    hold on;
    ee = 0.25e-01;
    plot([-1 1], [z0-ee z0-ee], '-k');
    hold on;
    plot([-1 1], [ze+ee ze+ee], '-k');

    Planes = getPlanes(input_multem.spec_atoms, input_multem.spec_lx, input_multem.spec_ly, inFP, input_multem.DimFP, input_multem.SeedFP);
    [nPlanes, ~] = size(Planes);
    for i = 1:nPlanes
        hold on;
        plot([-1 1], [Planes(i) Planes(i)], '-b');
        set(gca,'FontSize',10,'LineWidth',1,'PlotBoxAspectRatio',[0.75 1 1]);
    end
    [nAtoms, nslice, nPlanes]
end

