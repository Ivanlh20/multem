function [x] = ilm_rand(x_min, x_max, nr, nc)
    if(nargin<4)
        nc = 1;
    end
    if(nargin<3)
        nr = 1;
    end
        
    x = x_min+(x_max-x_min)*rand(nr, nc);
end