function[bb] = ilm_chk_bound(x, x_min, x_max)
    bb = (x_min<=x) && (x<x_max);
end