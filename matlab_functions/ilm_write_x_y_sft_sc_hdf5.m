function ilm_write_x_y_sft_sc_hdf5(fn, x, y, x_sft, x_sc, y_sft, y_sc)
    x_shape = size(x);
    y_shape = size(y);

    bb_y_sc = ~(nargin<7);
    bb_y_sft = ~(nargin<6);
    bb_x_sc = ~(nargin<5);
    bb_x_sft = ~(nargin<4);
 
    if bb_x_sft
        x_sft_shape = size(x_sft);
    end

    if bb_x_sc
        x_sc_shape = size(x_sc);
    end

    if bb_y_sft
        y_sft_shape = size(y_sft);
    end

    if bb_y_sc
        y_sc_shape = size(y_sc);
    end    
    
    bb = false;
    for ik=1:4
        try
            if exist(fn, 'file')
                delete(fn);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h5create(fn, '/x', x_shape, 'Datatype', class(x));
            h5create(fn, '/y', y_shape, 'Datatype', class(y));

            if bb_x_sft
                h5create(fn, '/x_sft', x_sft_shape, 'Datatype', class(x_sft));
            end

            if bb_x_sc
                h5create(fn, '/x_sc', x_sc_shape, 'Datatype', class(x_sc));
            end
            
            if bb_y_sft
                h5create(fn, '/y_sft', y_sft_shape, 'Datatype', class(y_sft));
            end

            if bb_y_sc
                h5create(fn, '/y_sc', y_sc_shape, 'Datatype', class(y_sc));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h5write(fn, '/x', x);
            h5write(fn, '/y', y);

            if bb_x_sft
                h5write(fn, '/x_sft', x_sft);
            end

            if bb_x_sc
                h5write(fn, '/x_sc', x_sc);
            end

            if bb_y_sft
                h5write(fn, '/y_sft', y_sft);
            end

            if bb_y_sc
                h5write(fn, '/y_sc', y_sc);
            end

            bb = true;
            break;
        catch
            pause(1);
        end
    end
    
    if ~bb
        disp(['data will not be saved: ', fn]); 
        if exist(fn, 'file')
            delete(fn);
        end
    end
end