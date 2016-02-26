function[]=show_crystal(fig, atoms)
    aZ = unique(atoms(:, 1)); 
    figure(fig); clf;
    for iZ=1:length(aZ)
        hold all;
        ii = find(atoms(:, 1)==aZ(iZ));
        plot3(atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'o');
    end;
    axis equal;
    view([1 1 1]);
end