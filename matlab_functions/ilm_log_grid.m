function[x] = ilm_log_grid(x_min, x_max, nx)
    dlnx = log(x_max/x_min)/(nx-1);
    x = x_min*exp((0:1:(nx-1))*dlnx).';
end