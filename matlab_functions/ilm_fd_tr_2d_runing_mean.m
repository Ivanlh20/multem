% find shifting between images
function [at_p] = ilm_fd_tr_2d_runing_mean(system_config, data, p, sigma_g, radius, bd_0, at_p0, n_it, b_fit, ik_ref)
    if(isstruct(data))
        disp('Error: The data must be as 3d stack')
        return;
    end
    
    [ny, nx, n_data] = ilm_size_data(data);
    
     if(nargin<10)
        ik_ref = 1; % It has to be implemented in the future
     end
     
     if(nargin<9)
        b_fit = true;
     end
    
    if(nargin<8)
        n_it = 1;
    end
    
    if(nargin<7)
        at_p0 = ilm_dflt_at_p(n_data);
    end
    
    if(nargin<6)
        bd_0 = zeros(1, 4);
    end

    % convert affine transformation to reference position
    at_p = at_p0;
    at_p(:,5:6) = at_p(:,5:6) - at_p(ik_ref, 5:6);
    
    % find shift between images
    dx = 1;
    dy = 1;
    for it = 1:n_it
        ik = 1;
        im_ik = double(data(:, :, ik));
        im_rm = ilc_tr_2d(im_ik, 1, 1, at_p(ik, 5:6), 0);
        tic;
        for ik = 2:n_data
            at_b = at_p(ik, :);
            bd_t = ilm_at_v_2_bd(nx, ny, bd_0, at_p0(ik-1, :));
            bd_ik = ilm_set_borders(at_b(5:6));
            bd = max(bd_t, bd_ik);
            
            [A, txy] = ilm_cvt_at_v_2_A_txy(at_b);
            im_ik = double(data(:, :, ik));
            at_p(ik, 5:6) = ilc_fd_tr_2d(system_config, im_rm/(ik-1), im_ik, dx, dy, A, txy, p, sigma_g, bd, 1, 0, b_fit, radius);

            im_rm = im_rm + ilc_tr_2d(im_ik, 1, 1, at_p(ik, 5:6), 0);

            if mod(ik, 10)==0
                t_c = toc;
                disp(['iter = ', num2str(it), ' calculating shift between images #', num2str(ik-1), '_', num2str(ik), '_time=', num2str(t_c, '%5.2f')]);
                tic;
            end
            
            if 0
               [i_min,i_max] = ilm_min_max(im_ik);
               figure(3);
               subplot(1, 3, 1);
               imagesc(im_rm/(ik-1), [i_min,i_max]);
               axis image off;
               colormap jet;
               title(['# =:', num2str(ik)])
               subplot(1, 3, 2);
               imagesc(im_ik,[i_min,i_max]);
               axis image off;
               colormap jet;
               title(['txy = [', num2str(at_p(ik, 5), 2), ', ', num2str(at_p(ik, 6), 2), ']'])
               subplot(1, 3, 3);
               plot(at_p(1:ik, 5), at_p(1:ik, 6), '-or')
            end
            % ilm_show_pcf(system_config, double(data(:, :, ik-1)), double(data(:, :, ik)), at_b, at_p(ik, :), p, sigma_g, bd);
        end
        
        % transform affine transformations related to the reference image   
        at_p(:,5:6) = at_p(:,5:6) - at_p(ik_ref, 5:6);
    end
end