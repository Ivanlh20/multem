function [x, y, x_sft, x_sc, y_sft, y_sc, y_dr] = ilm_read_x_y_sft_sc_y_dr_hdf5(fn)
    x = h5read(fn, '/x');
    y = h5read(fn, '/y');

    try
        x_sft = h5read(fn, '/x_sft');
    catch
        x_sft = 0.0;
    end

    try
        x_sc = h5read(fn, '/x_sc');
    catch
        x_sc = 1.0;
    end

    try
        y_sft = h5read(fn, '/y_sft');
    catch
        y_sft = 0.0;
    end

    try
        y_sc = h5read(fn, '/y_sc');
    catch
        y_sc = 1.0;
    end

    try
        y_dr = h5read(fn, '/y_dr');
    catch
        y_dr = 1.0;
    end   
end