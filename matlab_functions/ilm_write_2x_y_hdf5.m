function ilm_write_2x_y_hdf5(fn, x, x_1, y)
    x_shape = size(x);
    x_1_shape = size(x_1);
    y_shape = size(y);
    
    if exist(fn, 'file')
        try
            info = h5info(fn, '/x');
            bb = isequal(info.Dataspace.Size, x_shape);
            info = h5info(fn, '/x_1');
            bb = bb && isequal(info.Dataspace.Size, x_1_shape);
            info = h5info(fn, '/y');
            bb = bb && isequal(info.Dataspace.Size, y_shape);
            if ~bb
                delete(fn)
                
                h5create(fn, '/x', x_shape, 'Datatype', class(x))
                h5create(fn, '/x_1', x_1_shape, 'Datatype', class(x_1))
                h5create(fn, '/y', y_shape, 'Datatype', class(y))
            end
        catch
            delete(fn)
            
            h5create(fn, '/x', x_shape, 'Datatype', class(x))
            h5create(fn, '/x_1', x_1_shape, 'Datatype', class(x_1))
            h5create(fn, '/y', y_shape, 'Datatype', class(y))
        end
    else
        h5create(fn, '/x', x_shape, 'Datatype', class(x))
        h5create(fn, '/x_1', x_1_shape, 'Datatype', class(x_1))
        h5create(fn, '/y', y_shape, 'Datatype', class(y))
    end

    h5write(fn, '/x', x)
    h5write(fn, '/x_1', x_1)
    h5write(fn, '/y', y)
end