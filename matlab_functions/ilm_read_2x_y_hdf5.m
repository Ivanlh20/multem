function [x, x_1, y] = ilm_read_2x_y_hdf5(fn)
    x = h5read(fn, '/x');
    x_1 = h5read(fn, '/x_1');
    y = h5read(fn, '/y');
end