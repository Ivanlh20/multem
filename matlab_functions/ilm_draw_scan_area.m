% Draw area based p = [p_0;p_e]
function[] = ilm_draw_scan_area(fig, p, c)
    if(nargin<3)
        c = '-r';
    else
        c = ['-', erase(c, '-')];
    end
    lxy = diff(p, 1, 1);
    
    figure(fig);
    hold on;
    plot([p(1, 1), p(1, 1)+lxy(1)], [p(1, 2), p(1, 2)], c, ...
    [p(1, 1)+lxy(1), p(1, 1)+lxy(1)], [p(1, 2), p(1, 2)+lxy(2)], c, ...
    [p(1, 1)+lxy(1), p(1, 1)], [p(1, 2)+lxy(2), p(1, 2)+lxy(2)], c, ...
    [p(1, 1), p(1, 1)], [p(1, 2)+lxy(2), p(1, 2)], c);
    axis equal;
end