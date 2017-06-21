function [r] = ilm_linear_interp(x, p0, p1)
    y_min = min(p0(2), p1(2));
    y_max = max(p0(2), p1(2));
    m = 0;
    if(abs(p1(1)-p0(1))>1e-10)
        m = (p1(2)-p0(2))/(p1(1)-p0(1));
    end
    r = m*(x-p0(1))+p0(2);
    r = max(y_min, min(y_max, r));
end