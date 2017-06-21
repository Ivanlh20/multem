function [x] = ilm_randn(xm, xs, x_min, x_max)
    x = xm+xs*randn();
    while((x_min>x)||(x>x_max))
        x = xm+xs*randn();
    end
end