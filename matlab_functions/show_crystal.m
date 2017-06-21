function[]=show_crystal(fig, atoms)
    aZ = unique(atoms(:, 1)); 
    figure(fig); clf;
    for iZ=1:length(aZ)
        hold all;
        ii = find(atoms(:, 1)==aZ(iZ));
        plot3(atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'auto');
    end
    axis equal;
    xlim([min(atoms(:, 2)), max(atoms(:, 2))])
    ylim([min(atoms(:, 3)), max(atoms(:, 3))])
    zlim([min(atoms(:, 4))-2, max(atoms(:, 4))+2])
    view([1 1 1]);
    set(gca,'zdir','reverse');
end