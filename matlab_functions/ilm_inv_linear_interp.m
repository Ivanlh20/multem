function [r] = ilm_inv_linear_interp(x, p0, p1)
    x = 1/x;
    p0(1) = 1/p0(1);
    p1(1) = 1/p1(1);
    r = ilm_linear_interp(x, p0, p1);
end