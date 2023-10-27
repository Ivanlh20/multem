function [sz] = ilm_read_len_hdf5(fn, str_var, dim)
    if exist(fn, 'file')
        try
            info = h5info(fn, ['/', str_var]);
            sz = info.Dataspace.Size(dim);
        catch
            sz = 0; 
        end
    else
        sz = 0;
    end

end