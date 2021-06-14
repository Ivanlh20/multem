function [pp]= fcn_quad_spline(x, y)
    x = reshape(x, [], 1);
    y = reshape(y, [], 1);
    n = length(x);
    dx = diff(x);
    dy = diff(y);
    b = [0; 2*dy(1:(end-1))]./[1; dx(1:(end-1))];
    e = ones(n-1, 1);
    M = spdiags([e, e], [-1, 0], n-1, n-1);
    c_2 = M\b;
    c_1 = (dy - c_2.*dx)./dx.^2;
    c = [c_1, c_2, y(1:(n-1))];
    pp = mkpp(x, c);
end