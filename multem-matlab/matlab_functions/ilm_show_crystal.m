function[]=ilm_show_crystal(fig, atoms, bb_clf, bb_rz)
   if(nargin<4)
        bb_rz = true;
    end
    
    if(nargin<3)
        bb_clf = true;
    end
    
    aZ = unique(atoms(:, 1));
    
    if(fig>0)
        figure(fig);
    end
    
    if(bb_clf)
        clf;
    else
        hold on;
    end
    
    str = {};
    for iZ=1:length(aZ)
        hold all;
        ii = find(atoms(:, 1)==aZ(iZ));
        plot3(atoms(ii, 2), atoms(ii, 3), atoms(ii, 4), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'auto');
        str{iZ} = num2str(aZ(iZ));
    end
%  legend(str);
    axis equal;
    r_min = min(atoms(:, 2:4), [], 1);
    r_max = max(atoms(:, 2:4), [], 1);

    x_min = r_min(1)-1;
    x_max = r_max(1)+1;
    y_min = r_min(2)-1;
    y_max = r_max(2)+1;
    z_min = r_min(3)-1;
    z_max = r_max(3)+1;
    
    axis([x_min, x_max, y_min, y_max, z_min, z_max]);
    view([0 0 1]);
    
    if(bb_rz)
        set(gca, 'zdir', 'reverse');
    end
end