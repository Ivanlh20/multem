function[data, nr_grid, at_p] = ilm_gui_stem_reg(system_config, data_i, at_p0)
    global n_data nx ny x_c y_c sigma_g_max p fim
    global fz_t dx dy b_stop bb_nr p_rect_sel
    global f_xy fs_ax fs_ay fs_x_c fs_y_c
    global bb_shift_kpf bb_ctrl_kpf
    
    close all
    
    bb_shift_kpf = false;
    bb_ctrl_kpf = false;
    b_stop = false;
    
    if(nargin<3)
        n_data = ilm_nz_data(data_i);
        at_p0 = ilm_dflt_at_p(n_data);
    end

    fy = 0.60;
    fx = 0.75;
    dm = get(0, 'MonitorPositions');
    dm = dm(1, :);
    w = fx*dm(3);
    h = fy*dm(4);
    x0 = (1-fx)*w/2;
    y0 = (1-fy)*h/2;

    fh = figure(1); clf;
    set(fh, {'Name', 'Visible', 'units', 'position', 'CloseRequestFcn'}, {'STEM registration -- Software created by Ivan Lobato: Ivanlh20@gmail.com', 'on', 'pixels', [x0, y0, x0+w, y0+h], @cb_close})
    set(fh, 'MenuBar', 'none');
    set(fh, 'ToolBar', 'none');    
    movegui(fh, 'center');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ui_pnl_opt = uipanel('Parent', fh, 'Title', 'Options', 'Units', 'normalized', ..., 
    'FontSize', 12, 'BackgroundColor', 'white', 'Position', [0.1 0.725 0.70 0.26]);

    pnl_x_0 = 12;
    pnl_y = 3 + (4:-1:0)*35;
    pnl_d = 4;
    pnl_h = 24;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1, 2,  3,   4,  5, 6,  7,   8,  9, 10, 11, 12, 13, 14
    pnl_w = [80, 70, 70, 50, 50, 50, 65, 40, 40, 45, 55, 55, 75, 45];
    pnl_n = length(pnl_w);
    pnl_x = pnl_x_0*ones(1, pnl_n);
    
    for it=2:pnl_n
        pnl_x(it) = pnl_x(it-1) + pnl_d + pnl_w(it-1);
    end

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'Rigid Reg.', 'Position', [pnl_x(1) pnl_y(1) pnl_w(1) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Reset txy', 'Position', [pnl_x(2) pnl_y(1) pnl_w(2) pnl_h], 'Callback', @cb_reset_txy);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Reset A', 'Enable', 'off', 'Position', [pnl_x(3) pnl_y(1) pnl_w(3) pnl_h], 'Callback', @cb_reset_A);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'Crt. tx', 'Position', [pnl_x(4) pnl_y(1) pnl_w(4) pnl_h], 'Callback', @(src, ev)cb_crt_txy(1));

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'Crt. ty', 'Position', [pnl_x(5) pnl_y(1) pnl_w(5) pnl_h], 'Callback', @(src, ev)cb_crt_txy(2));

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'Crt. txy', 'Position', [pnl_x(6) pnl_y(1) pnl_w(6) pnl_h], 'Callback', @(src, ev)cb_crt_txy(3));

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'Crt. A_txy', 'Enable', 'off', 'Position', [pnl_x(7) pnl_y(1) pnl_w(7) pnl_h], 'Callback', @cb_crt_A_txy);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'n_it', 'Position', [pnl_x(8) pnl_y(1) pnl_w(8) pnl_h]);
    
    ui_edt_rg_n_it = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '1', 'Position', [pnl_x(9) pnl_y(1) pnl_w(9) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'idx_ref', 'Position', [pnl_x(10) pnl_y(1) pnl_w(10) pnl_h]);
    
    ui_pu_idx_ref = uicontrol('Parent', ui_pnl_opt, 'Style', 'popup', ...
    'String', num2cell(1:n_data), 'Value', 1, ...
    'Position', [pnl_x(11) pnl_y(1) pnl_w(11) pnl_h], 'Callback', @cb_idx_ref);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'Preproc', 'Position', [pnl_x(12) pnl_y(1) pnl_w(12) pnl_h]);
    
    ui_pu_preproc = uicontrol('Parent', ui_pnl_opt, 'Style', 'popup', ...
    'String', {'Raw', 'Laplacian', 'LCWT'}, 'Value', 2, ...
    'Position', [pnl_x(13) pnl_y(1) pnl_w(13) pnl_h], 'Callback', @cb_preproc);

    ui_edt_preproc_sigma = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '4.0', 'Enable', 'off', 'Position', [pnl_x(14) pnl_y(1) pnl_w(14) pnl_h], 'Callback', @cb_show_sel_ui_lb_data);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1, 2,  3,  4,  5,   6,  7,  8,  9, 10, 11
    pnl_w = [80, 80, 80, 45, 55, 45, 55, 45, 55, 45, 60];
    pnl_n = length(pnl_w);
    pnl_x = pnl_x_0*ones(1, pnl_n);

    for it=2:pnl_n
        pnl_x(it) = pnl_x(it-1) + pnl_d + pnl_w(it-1);
    end  

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'Non-Rigid Reg.', 'Position', [pnl_x(1) pnl_y(2) pnl_w(1) pnl_h]);
    
    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Assign A-txy', 'Position', [pnl_x(2) pnl_y(2) pnl_w(2) pnl_h], 'Callback', @cb_assign_A_txy);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'Crt. nr_grid', 'Position', [pnl_x(3) pnl_y(2) pnl_w(3) pnl_h], 'Callback', @cb_crt_nr_grid);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'n_it_i', 'Position', [pnl_x(4) pnl_y(2) pnl_w(4) pnl_h]);
    
    ui_edt_nrg_n_it_i = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '10', 'Position', [pnl_x(5) pnl_y(2) pnl_w(5) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'n_it_o', 'Position', [pnl_x(6) pnl_y(2) pnl_w(6) pnl_h]);
    
    ui_edt_nrg_n_it_o = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '8', 'Position', [pnl_x(7) pnl_y(2) pnl_w(7) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'alpha', 'Position', [pnl_x(8) pnl_y(2) pnl_w(8) pnl_h]);
    
    ui_edt_alpha = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '0.5', 'Position', [pnl_x(9) pnl_y(2) pnl_w(9) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'Opt.', 'Position', [pnl_x(10) pnl_y(2) pnl_w(10) pnl_h]);
    
    ui_pu_opt = uicontrol('Parent', ui_pnl_opt, 'Style', 'popup', ... 
    'String', {'Mean', 'Poly 1', 'Poly 2', 'orig.'}, 'Position', [pnl_x(11) pnl_y(2) pnl_w(11) pnl_h]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1,   2,  3,  4,  5,  6,   7
    pnl_w = [270, 50, 50, 50, 65, 100, 110];
    pnl_n = length(pnl_w);
    pnl_x = pnl_x_0*ones(1, pnl_n);   

    for it=2:pnl_n
        pnl_x(it) = pnl_x(it-1) + pnl_d + pnl_w(it-1);
    end  
    
    ui_bg_opt = uibuttongroup('Parent', ui_pnl_opt, 'units', 'Pixels', ...
        'Position', [pnl_x(1) pnl_y(3) pnl_w(1) pnl_h], 'SelectionChangedFcn', @cb_execute_selection);

    ui_rb_select = uicontrol('Parent', ui_bg_opt, 'Style', 'radiobutton', ...
        'String', 'select', 'Position', [pnl_x_0 0 50 pnl_h]);

    ui_rb_add_bd = uicontrol('Parent', ui_bg_opt, 'Style', 'radiobutton', ...
        'String', 'add bd', 'Position', [pnl_x_0+(60+pnl_d) 0 60 pnl_h]);

    ui_rb_resize = uicontrol('Parent', ui_bg_opt, 'Style', 'radiobutton', ...
        'String', 'resize', 'Position', [pnl_x_0+(120+2*pnl_d) 0 60 pnl_h]);
    
    ui_rb_crop = uicontrol('Parent', ui_bg_opt, 'Style', 'radiobutton', ...
        'String', 'crop', 'Position', [pnl_x_0+(180+3*pnl_d) 0 60 pnl_h]);

    ui_edt_sampling_x = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '0', 'Position', [pnl_x(2) pnl_y(3) pnl_w(2) pnl_h]);

    ui_edt_sampling_y = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '0', 'Position', [pnl_x(3) pnl_y(3) pnl_w(3) pnl_h]);

    ui_pb_execute = uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Exec.', 'Position', [pnl_x(4) pnl_y(3) pnl_w(4) pnl_h], 'Callback', @cb_execute);

    ui_cb_reg_typ = uicontrol('Parent', ui_pnl_opt, 'Style', 'checkbox', ... 
    'String', 'Rig. Reg.', 'Value', 1, 'Position', [pnl_x(5) pnl_y(3) pnl_w(5) pnl_h], 'Callback', @cb_reg_typ);

    ui_pb_show_data = uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Show rg data', 'Position', [pnl_x(6) pnl_y(3) pnl_w(6) pnl_h], 'Callback', @cb_show_data);

    ui_pb_show_avg_data = uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'show avg rg data', 'Position', [pnl_x(7) pnl_y(3) pnl_w(7) pnl_h], 'Callback', @cb_show_avg_data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % 1,   2,  3,  4,  5,  6,  7,  8, 9,  10,  11
    pnl_w = [75, 75, 60, 70, 60, 60, 60, 65, 165, 60, 60];
    pnl_n = length(pnl_w);
    pnl_x = pnl_x_0*ones(1, pnl_n);   

    for it=2:pnl_n
        pnl_x(it) = pnl_x(it-1) + pnl_d + pnl_w(it-1);
    end 
    
    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'Reset data', 'Position', [pnl_x(1) pnl_y(4) pnl_w(1) pnl_h], 'Callback', @cb_reset_data);
    
    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ...
    'String', 'Sigma_g(px.)', 'Position', [pnl_x(2) pnl_y(4) pnl_w(2) pnl_h]);
    
    ui_edt_sigma = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ...
    'String', num2str(sigma_g_max/2, '%5.2f'), 'Position', [pnl_x(3) pnl_y(4) pnl_w(3) pnl_h]);

    ui_cb_rcfft = uicontrol('Parent', ui_pnl_opt, 'Style', 'checkbox', ... 
    'String', 'Recal. fft', 'Value', 1, 'Position', [pnl_x(4) pnl_y(4) pnl_w(4) pnl_h]);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ...
    'String', 'Exp. FFT', 'Position', [pnl_x(5) pnl_y(4) pnl_w(5) pnl_h]);
    
    ui_pu_exp = uicontrol('Parent', ui_pnl_opt, 'Style', 'popup', ...
    'String', {'0.100', '0.125', '0.20', '0.25', '0.5', '0.75'}, 'Value', 3, ...
    'Position', [pnl_x(6) pnl_y(4) pnl_w(6) pnl_h], 'Callback', @(src, ev)fcn_show_ax_f_1(false));

    uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ... 
    'String', 'Time(s.)', 'Position', [pnl_x(7) pnl_y(4) pnl_w(7) pnl_h]);
    
    ui_edt_time = uicontrol('Parent', ui_pnl_opt, 'Style', 'edit', ... 
    'String', '0.05', 'Position', [pnl_x(8) pnl_y(4) pnl_w(8) pnl_h]);

    ui_txt_msg = uicontrol('Parent', ui_pnl_opt, 'Style', 'text', ...
    'String', '00.0%', 'Position', [pnl_x(9) pnl_y(4) pnl_w(9) pnl_h], 'FontSize', 13, 'FontWeight', 'bold');

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ...
    'String', 'stop', 'Position', [pnl_x(10) pnl_y(4) pnl_w(10) pnl_h], 'Callback', @cb_stop);

    uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'String', 'close', 'Position', [pnl_x(11) pnl_y(4) pnl_w(11) pnl_h], 'Callback', @cb_close);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pnl_w = [850];
    pnl_n = length(pnl_w);
    pnl_x = pnl_x_0*ones(1, pnl_n);
    
    for it=2:pnl_n
        pnl_x(it) = pnl_x(it-1) + pnl_d + pnl_w(it-1);
    end 
    
    ui_pb_data_info = uicontrol('Parent', ui_pnl_opt, 'Style', 'pushbutton', ... 
    'Enable', 'inactive', 'String', ' ', 'Position', [pnl_x(1) pnl_y(5) pnl_w(1) pnl_h]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d_f = 0.010;
    h_f = 0.55;
    w_1 = 0.16;
    w_2 = h_f*h/w;
    w_3 = w_2;
    w_4 = w_3;
  
    x_1 = (1-(w_1+w_2+w_3+w_4+3*d_f))/2;
    x_2 = x_1+w_1+d_f;
    x_3 = x_2+w_2+d_f;
    x_4 = x_3+w_3+d_f;
    y_0 = 0.10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ui_lb_data = uicontrol('Style', 'listbox', 'Units', 'normalized', ...
        'Position', [x_1 y_0 w_1 h_f], 'Callback', @cb_show_sel_ui_lb_data);
    
    ax_f(1) = axes('Position', [x_2 y_0 w_2 h_f], 'Visible', 'off');
    ax_f(2) = axes('Position', [x_3 y_0 w_3 h_f], 'Visible', 'off');
    ax_f(3) = axes('Position', [x_4 y_0 w_4 h_f], 'Visible', 'off');
    
    fcn_init(system_config, data_i, at_p0);
    
    fh.WindowButtonMotionFcn = @cb_mouse_move;
    fh.WindowButtonDownFcn = @cb_mouse_click;
    fh.WindowScrollWheelFcn = @cb_zoom;
    fh.WindowKeyPressFcn = @cb_press_key;
    fh.WindowKeyReleaseFcn = @cb_release_key;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uiwait(fh)
    
    function [typ] = fcn_get_activate_typ()
        typ = ilm_ifelse(ui_rb_add_bd.Value, 1, 0);
        typ = ilm_ifelse(ui_rb_resize.Value, 2, typ);
        typ = ilm_ifelse(ui_rb_crop.Value, 3, typ);
   end
        
    function fcn_activate_opt(typ)
        if(typ==0) % select
            ui_edt_sampling_x.Enable = 'off';
            ui_edt_sampling_y.Enable = 'off';
            
            ui_pb_execute.Enable = 'off';
        elseif(typ==1) % add bd
            ui_edt_sampling_x.String = num2str(0);
            ui_edt_sampling_x.Enable = 'on';
            
            ui_edt_sampling_y.String = num2str(0);
            ui_edt_sampling_y.Enable = 'on';
            
            ui_pb_execute.Enable = 'on';
            ui_pb_execute.Enable = 'on';
        elseif(typ==2) % resize 
            if(str2double(ui_edt_sampling_x.String)<64)
                ui_edt_sampling_x.String = num2str(nx);
            end
            ui_edt_sampling_x.Enable = 'on';
            
            if(str2double(ui_edt_sampling_y.String)<64)
                ui_edt_sampling_y.String = num2str(ny);
            end
            ui_edt_sampling_y.Enable = 'off';
            
            ui_pb_execute.Enable = 'on';
        elseif(typ==3) % crop
            ui_edt_sampling_x.Enable = 'off';
            ui_edt_sampling_y.Enable = 'off';
            
            ui_pb_execute.Enable = 'on';
        end
    end

    function cb_execute_selection(source, event)
        typ = fcn_get_activate_typ();
        fcn_activate_opt(typ)
    end

    function cb_execute(source, event)
        fcn_msn_start();
        
        typ = fcn_get_activate_typ();
        
        if(typ==1) % add bd
            bd_x = round(str2double(ui_edt_sampling_x.String));
            bd_x = ilm_pn_border(nx, bd_x, 1);
            nx_r = nx+2*bd_x;
            ix_0 = bd_x + 1;
            ix_e = ix_0 + nx-1;
            ax = ix_0:ix_e;
            
            bd_y = round(str2double(ui_edt_sampling_y.String));
            bd_y = ilm_pn_border(ny, bd_y, 1);
            ny_r = ny + 2*bd_y;
            iy_0 = bd_y + 1;
            iy_e = iy_0 + ny-1;
            ay = iy_0:iy_e;
            
            data_t = zeros(ny_r, nx_r, n_data);
            % data_t(ay, ax, :) = data;
            for ik=1:n_data
                im_ik = double(data(:, :, ik));
                data_t(:, :, ik) = mean(im_ik(:));
                data_t(ay, ax, ik) = im_ik;

                ui_txt_msg.String = [num2str(100*ik/n_data, '%4.1f'), ' %'];
                pause(0.001);
            end
            
            p_rect_sel_0 = p_rect_sel + [bd_x, bd_y];
            
            fcn_init(system_config, data_t, at_p, p_rect_sel_0);
            
            ui_edt_sampling_x.String = num2str(0);
            ui_edt_sampling_y.String = num2str(0);

        elseif(typ==2) % resize 
            nx_r = round(str2double(ui_edt_sampling_x.String));
            
            if(nx_r==0)
                return;
            end
            ny_r = round(nx_r*ny/nx);
            nx_r = ilm_pn_fact(nx_r);

            data_t = zeros(ny_r, nx_r, n_data);
            for ik=1:n_data
                im_ik = double(data(:, :, ik));
                data_t(:, :, ik) = max(0, imresize(im_ik, [ny_r, nx_r], 'bicubic'));

                ui_txt_msg.String = [num2str(100*ik/n_data, '%4.1f'), ' %'];
                pause(0.001);
            end

            at_p(:, 5:6) = at_p(:, 5:6).*[nx_r/nx, ny_r/ny];
            if(mod(ny_r, 2)==1)
                data_t = data_t(1:(end-1), :, :);
            end
            
            fcn_init(system_config, data_t, at_p);
            
            ui_edt_sampling_x.String = num2str(0);
            ui_edt_sampling_y.String = num2str(0);
        elseif(typ==3) % crop
            ax = p_rect_sel(1, 1):p_rect_sel(2, 1);
            ay = p_rect_sel(1, 2):p_rect_sel(2, 2);

            for ik=1:n_data
                [A, txy] = ilm_cvt_at_v_2_A_txy(at_p(ik, :));
                im_ik = double(data(:, :, ik));
                data(:, :, ik) = ilc_tr_at_2d(system_config, im_ik, dx, dy, A, txy, 3, 0);
                ui_txt_msg.String = [num2str(100*ik/n_data, '%4.1f'), ' %'];
                pause(0.001);
            end
            
            at_p = ilm_dflt_at_p(n_data);
            
            fcn_init(system_config, data(ay, ax, :), at_p);
        end
        
        fcn_msn_end();
    end

    function s = fcn_str_data_info(p_rect_sel)
        s = ['[nx, ny] = [', num2str(nx), ', ', num2str(ny), ']'];
        s = [s, ' - [ix_0, ix_e] = [', num2str(p_rect_sel(1, 1)), ', ', num2str(p_rect_sel(2, 1)), ']'];
        s = [s, ' - [iy_0, iy_e] = [', num2str(p_rect_sel(1, 2)), ', ', num2str(p_rect_sel(2, 2)), ']'];
        s = [s, ' - [nx_s, ny_s] = [', num2str(abs(p_rect_sel(2, 1)-p_rect_sel(1, 1))+1), ', ', num2str(abs(p_rect_sel(2, 2)-p_rect_sel(1, 2))+1), ']'];
    end

    function fcn_init(system_config, data_i, at_p_0, p_rect_sel_0)  
        fcn_msn_start();

        bb_nr = false;
        
        data = data_i;
        [ny, nx, n_data] = ilm_size_data(data);
        
        if(nx>ny)
            f_xy = [1, nx/ny];
        else
            f_xy = [ny/nx, 1];
        end
        
        fs_ax = (1:nx)*f_xy(1);
        fs_ay = (1:ny)*f_xy(2);
        fs_x_c = nx*f_xy(1)/2 + 1;
        fs_y_c = ny*f_xy(2)/2 + 1;
        
        dx = 1;
        dy = 1;
        
        x_c = nx/2 + 1;
        y_c = ny/2 + 1;
        
        if(nargin<4)
            p_rect_sel_0 = [1, 1; nx, ny];
        end
        
        p_rect_sel = p_rect_sel_0;
        
        at_p = at_p_0;

        % set sigma max
        sigma_g_max = min(nx/2, ny/2)-1;

        % set initial sigma
        fim = fcn_mean_abs_fdata_at(system_config, data, at_p);
        sigma_gt = ilc_info_lim_2d(fim);
        
        ui_edt_sampling_x.String = num2str(nx);
        ui_edt_sampling_y.String = num2str(ny);
        ui_edt_sigma.String = num2str(sigma_gt, '%5.2f');
        ui_pb_data_info.String = fcn_str_data_info(p_rect_sel);
        
        % set initial zoom    
        fz_t(1) = sigma_g_max/(2*sigma_gt);
        fz_t(2) = 1;
        fz_t(3) = sigma_g_max/min(sigma_g_max/2, max(sigma_g_max/5, 6*ilm_sigma_g_2_sigma_r(nx, ny, sigma_gt)));

        % plot fim
        fcn_show_ax_f_1(false);
        
        axes(ax_f(1));
        ax = gca;
        ax.Toolbar.Visible = 'off';
        
        % set initial lb_data values
        fcn_set_at_p_2_ui_lb_data(at_p);
        
        axes(ax_f(2));
        ax = gca;
        ax.Toolbar.Visible = 'off';
        
        axes(ax_f(3));
        ax = gca;
        ax.Toolbar.Visible = 'off';
        
        ui_lb_data.Value = 1;
        cb_show_sel_ui_lb_data();
        
        fcn_msn_end();
    end

    function [image_av] = fcn_mean_abs_fdata_at(system_config, data, at_p)
        [ny, nx, n_data] = ilm_size_data(data);

        n_idx = max(2, round(n_data/2));
        a_idx = randperm(n_data, n_idx);

        image_av = zeros(ny, nx);
        radius = 0.85*min(nx, ny)/2;
        fbw = ilc_wd_butwth([nx, 1], [1, 1], 8, radius);
        fbw = fbw.*ilc_wd_butwth([1, ny], [1, 1], 8, radius);

        for ik = a_idx
            [A, txy] = ilm_cvt_at_v_2_A_txy(at_p(ik, :));
            image_ik = double(data(:, :, ik));
            [image_ik, bg_ik] = ilc_tr_at_2d(system_config, image_ik, dx, dy, A, txy, 3, 0);
            image_ik = (image_ik-0.999*bg_ik).*fbw;

            image_av = image_av + abs(fft2(image_ik));
        end

        image_av = ifftshift(image_av/n_idx);
    end

    function fcn_msn_start()
        ui_txt_msg.String = 'Processing';
        drawnow;
    end
 
    function fcn_msn_end()
        ui_txt_msg.String = 'Done';
        drawnow;
    end

    function fcn_set_at_p_2_ui_lb_data(at_p)
        n_at_p = size(at_p, 1);
        items{1, n_at_p} = [];
        UserData{1, n_at_p} = [];
        for ik=1:n_at_p
            s_at_p_ik = ['[', num2str(at_p(ik, 3), '%4.3f'), ', ', num2str(at_p(ik, 4), '%4.3f'), ...
                ', ', num2str(at_p(ik, 5), '%5.2f'), ', ', num2str(at_p(ik, 6), '%5.2f'), ']'];

            items{ik} = ['ik = ', num2str(ik, '%.3d'), ', at_p = ' s_at_p_ik];
            UserData{ik} = ik;
        end
        ui_lb_data.String = items;
        ui_lb_data.UserData = UserData;
    end

    function[im] = fcn_rescale_img(im_i)
        alpha = str2double(ui_pu_exp.String{ui_pu_exp.Value});
        if(~isnumeric(alpha))
            alpha = 0.20;
        end
        
        alpha = max(0.01, abs(alpha));
        
        im = abs(im_i).^alpha;
    end

    function fcn_set_nr_grid_0(n_data, nx, ny)
        nr_grid.x = zeros(ny, nx, n_data, 'single');
        nr_grid.y = zeros(ny, nx, n_data, 'single');
    end

    function st_opt = fcn_st_preproc_st_opt()
        sigma = str2double(ui_edt_preproc_sigma.String);
        if(~isnumeric(sigma))
            sigma = 4.0;
        end
        st_opt.preproc_opt = ui_pu_preproc.Value;
        st_opt.preproc_sigma =  max(1.0, abs(sigma));  
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fcn_show_ax_f_1(b_fimage_rc)
        if(nargin==0)
            b_fimage_rc = true;
        end
        
        if(b_fimage_rc)
            fim = fcn_mean_abs_fdata_at(system_config, data, at_p);
        end
        
        % set axes
        axes(ax_f(1));
        
        % plot fimage
        imagesc(fs_ax, fs_ay, fcn_rescale_img(fim));
        axis image off;
        colormap bone;

        % plot circle
        sigma_gt = str2double(ui_edt_sigma.String);
        ilm_plot_circle(fs_x_c, fs_y_c, sigma_gt, 'red', 'f');

        zoom(fz_t(1));
    end

    function fcn_show_ax_f_2(im, title_str)
        % set axes
        axes(ax_f(2));
        imagesc(im);
        axis image off;
        colormap bone;
        title(title_str);
        
        bd = ilm_calc_borders_using_at_v(at_p);
        p_rect_sel_atp = fcn_bd_2_p_rect_sel(bd);
        ilm_plot_rectangle(p_rect_sel_atp, 'blue', 'f');
        
        ilm_plot_rectangle(p_rect_sel, 'red', 'f');

        zoom(fz_t(2));
    end

    function fcn_show_ax_f_3(im, title_str)
        % set axes
        axes(ax_f(3));
        imagesc(im);
        axis image off;
        colormap bone;
        title(title_str);

        zoom(fz_t(3));
    end

    function fcn_plot_data_rg(im_r, A_r, txy_r, im_s, A_s, txy_s, p, sigma_gt, bd, title_ax_2, title_ax_3)
        st_opt = fcn_st_preproc_st_opt();
        
        [im_rp, im_sp] = ilm_pcf_data_tf(im_r, im_s, st_opt.preproc_opt, st_opt.preproc_sigma);
        
        im_rp = ilc_tr_at_2d(system_config, im_rp, dx, dy, A_r, txy_r, 3, 0);

        pcf = ilc_pcf_2d(system_config, im_rp, im_sp, dx, dy, A_s, txy_s, p, sigma_gt, bd, 3, 0);
        
        im_r_plot = ilc_tr_at_2d(system_config, im_r, dx, dy, A_r, txy_r, 3, 0);
        fcn_show_ax_f_2(im_r_plot, title_ax_2);

        fcn_show_ax_f_3(pcf, title_ax_3);
    end
 
    function fcn_plot_data_nrg(im_r, im_s, p, sigma_gt, bd, title_ax_2, title_ax_3)
       
        A = [1 0; 0 1];
        txy = [0; 0];

        st_opt = fcn_st_preproc_st_opt();
        [im_rp, im_sp] = ilm_pcf_data_tf(im_r, im_s, st_opt.preproc_opt, st_opt.preproc_sigma);

        pcf = ilc_pcf_2d(system_config, im_rp, im_sp, dx, dy, A, txy, p, sigma_gt, bd, 3, 0);

        fcn_show_ax_f_2(im_r, title_ax_2);

        fcn_show_ax_f_3(pcf, title_ax_3);
    end

    function cb_show_sel_ui_lb_data(source, event)
        sigma_gt = str2double(ui_edt_sigma.String);
        
        if(sigma_gt>sigma_g_max)
            return
        end
        
        p = 0.05;
        bd = fcn_p_rect_sel_2_bd(p_rect_sel);
        
        fcn_msn_start();
        
        ik_r = ui_lb_data.UserData{ui_lb_data.Value};
        ik_s = ilm_ifelse(ik_r==n_data, ik_r-1, ik_r+1);
        
        title_ax_2 = ['Image # ', num2str(ik_r)];
        title_ax_3 = ['Pcf #= ', num2str(ik_r), ' - ', num2str(ik_s)];

        if ui_cb_reg_typ.Value
            im_r = double(data(:, :, ik_r));
            [A_r, txy_r] = ilm_cvt_at_v_2_A_txy(at_p(ik_r, :));

            im_s = double(data(:, :, ik_s));
            [A_s, txy_s] = ilm_cvt_at_v_2_A_txy(at_p(ik_s, :));

            fcn_plot_data_rg(im_r, A_r, txy_r, im_s, A_s, txy_s, p, sigma_gt, bd, title_ax_2, title_ax_3);
        else
            [Rx_i, Ry_i] = meshgrid(0:(nx-1), 0:(ny-1));
            
            Rx = Rx_i + double(nr_grid.x(:, :, ik_r));
            Ry = Ry_i + double(nr_grid.y(:, :, ik_r));
            im_r = ilc_intrpl_rn_2d(system_config, double(data(:, :, ik_r)), Rx, Ry);
            
            Rx = Rx_i + double(nr_grid.x(:, :, ik_s));
            Ry = Ry_i + double(nr_grid.y(:, :, ik_s));
            im_s = ilc_intrpl_rn_2d(system_config, double(data(:, :, ik_s)), Rx, Ry);
            
            fcn_plot_data_nrg(im_r, im_s, p, sigma_gt, bd, title_ax_2, title_ax_3)
        end
        
        fcn_msn_end();
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fcn_show_rg_data(sigma_gt)
        p = 0.05;
        bd = fcn_p_rect_sel_2_bd(p_rect_sel);
        for ik=2:n_data
            if(b_stop)
                break;
            end
            
            title_ax_2 = ['Image # ', num2str(ik-1)];
            title_ax_3 = ['Pcf #= ', num2str(ilm_set_bound(ik-1, 1, n_data+1)), ' - ', num2str(ik)];

            im_r = double(data(:, :, ik-1));
            [A_r, txy_r] = ilm_cvt_at_v_2_A_txy(at_p(ik-1, :));
        
            im_s = double(data(:, :, ik));
            [A_s, txy_s] = ilm_cvt_at_v_2_A_txy(at_p(ik, :));
            
            fcn_plot_data_rg(im_r, A_r, txy_r, im_s, A_s, txy_s, p, sigma_gt, bd, title_ax_2, title_ax_3)

            ui_txt_msg.String = [num2str(100*ik/n_data, '%4.1f'), ' %'];

            tm = str2double(ui_edt_time.String);
            tm = ilm_ifelse(isnan(tm), 0.075, tm);
            pause(tm);
        end
        b_stop = false;
    end

    function fcn_show_nr_data(sigma_gt)
        if(~bb_nr)
            nr_grid = ilm_at_p_2_nr_grid(at_p, nx, ny);
        end
        
        p = 0.05;
        bd = fcn_p_rect_sel_2_bd(p_rect_sel);
        
        [Rx_i, Ry_i] = meshgrid(0:(nx-1), 0:(ny-1));
        
        for ik=2:n_data   
            if(b_stop)
                break;
            end
            
            title_ax_2 = ['Image # ', num2str(ik-1)];
            title_ax_3 = ['Pcf #= ', num2str(ilm_set_bound(ik-1, 1, n_data+1)), ' - ', num2str(ik)];
            
            Rx = Rx_i + double(nr_grid.x(:, :, ik-1));
            Ry = Ry_i + double(nr_grid.y(:, :, ik-1));
            im_r = ilc_intrpl_rn_2d(system_config, double(data(:, :, ik-1)), Rx, Ry);
            
            Rx = Rx_i + double(nr_grid.x(:, :, ik));
            Ry = Ry_i + double(nr_grid.y(:, :, ik));
            im_s = ilc_intrpl_rn_2d(system_config, double(data(:, :, ik)), Rx, Ry);
            
            fcn_plot_data_nrg(im_r, im_s, p, sigma_gt, bd, title_ax_2, title_ax_3)

            ui_txt_msg.String = [num2str(100*ik/n_data, '%4.1f'), ' %'];

            tm = str2double(ui_edt_time.String);
            tm = ilm_ifelse(isnan(tm), 0.075, tm);
            pause(tm);
        end
        b_stop = false;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cb_reset_txy(source, event)
        at_p(:, 5:6) = at_p0(:, 5:6);

        % fill in ui_lb_data
        fcn_set_at_p_2_ui_lb_data(at_p);
            
        % select item in ui_lb_data
        cb_show_sel_ui_lb_data();

        % plot fimage
        fcn_show_ax_f_1(true);
    end

    function cb_reset_A(source, event)
        at_p(:, 1) = 1;
        at_p(:, 2:3) = 0;
        at_p(:, 4) = 1;
         
        % fill in ui_lb_data
        fcn_set_at_p_2_ui_lb_data(at_p);
            
        % select item in ui_lb_data
        cb_show_sel_ui_lb_data();

        % plot fimage
        fcn_show_ax_f_1(true);
    end

    function cb_assign_A_txy(source, event)
        fcn_msn_start();
        
        nr_grid = ilm_at_p_2_nr_grid(at_p, nx, ny);
        cb_show_sel_ui_lb_data();
        
        fcn_msn_end();
        
        bb_nr = true;
    end

    function cb_reset_data(source, event)
        data = data_i;
        n_data = ilm_nz_data(data);
        [ny, nx, n_data] = ilm_size_data(data);
        at_p = at_p0;

        % reset nr_grid
        fcn_set_nr_grid_0(n_data, nx, ny);
    
        % fill in ui_lb_data
        fcn_set_at_p_2_ui_lb_data(at_p);
 
        % plot fimage
        fcn_show_ax_f_1(true);

        % select item in ui_lb_data
        cb_show_sel_ui_lb_data();

        p_rect_sel = [1, 1; nx, ny];

        bb_nr = false;
    end

    function cb_stop(source, event)
        b_stop = true;
    end

    function cb_preproc(source, event)
        ui_edt_preproc_sigma.Enable = ilm_ifelse(ui_pu_preproc.Value==3, 'on', 'off');
        cb_show_sel_ui_lb_data(source, event);
    end

    function cb_idx_ref(source, event)
        idx_ref = ui_pu_idx_ref.Value;
        at_p(:, 5:6) = at_p(:, 5:6) - at_p(idx_ref, 5:6);
        cb_show_sel_ui_lb_data(source, event);
    end

    function cb_reg_typ(source, event)
        if source.Value
            ui_pb_show_data.String = 'Show rg data';
            ui_pb_show_avg_data.String = 'Show avg rg data';
        else
            ui_pb_show_data.String = 'Show nrg data';
            ui_pb_show_avg_data.String = 'Show avg nrg data';
        end

        cb_show_sel_ui_lb_data();
    end

    function cb_show_data(source, event)
        sigma_gt = str2double(ui_edt_sigma.String);

        if(sigma_gt>sigma_g_max)
            return
        end

        fcn_msn_start();
        
        if ui_cb_reg_typ.Value
            fcn_show_rg_data(sigma_gt);
        else
            fcn_show_nr_data(sigma_gt);
        end
        
        cb_show_sel_ui_lb_data();
        
        fcn_msn_end();
    end

    function cb_show_avg_data(source, event)    
        fcn_msn_start();
        
        if ui_cb_reg_typ.Value
            image_avg = ilm_mean_data_at(system_config, data, at_p, 3);
            fig_n = 2;
            title_str = 'Average image rigid registration';
        else
            image_avg = ilm_mean_data_nr_grid(system_config, data, nr_grid, 3);
            fig_n = 3;
            title_str = 'Average image non-rigid registration';
        end
        
        % set constant borders
        
        bd = ilm_calc_borders_using_at_v(at_p);
        image_avg = ilc_repl_bdr(image_avg, bd, 3);
    
        figure(fig_n);
        imagesc(image_avg);
        axis image off;
        colormap bone;
        title(title_str);
        
        fcn_msn_end();
        
        figure(fig_n);
    end

    function bd = fcn_p_rect_sel_2_bd(p_rect_sel)
        bd = [p_rect_sel(1, 1)-1, nx-p_rect_sel(2, 1), p_rect_sel(1, 2)-1, ny-p_rect_sel(2, 2)];
        bd(1) = max(0, bd(1));
        bd(2) = min(nx, bd(2));
        bd(3) = max(0, bd(3));
        bd(4) = min(ny, bd(4));
    end

    function p_rect = fcn_bd_2_p_rect_sel(bd)
        p_rect = [bd(1)+1, nx-bd(2); bd(3)+1, ny-bd(4)].';
        p_rect(1, 1) = max(1, p_rect(1, 1));
        p_rect(2, 1) = min(nx, p_rect(2, 1));
        p_rect(1, 2) = max(1, p_rect(1, 2));
        p_rect(2, 2) = min(ny, p_rect(2, 2));
    end

    function bd = fcn_atp_p_rect_sel_2_bd(at_p, p_rect_sel)
        bd = ilm_calc_borders_using_at_v(at_p);
        bd_user = [p_rect_sel(1, 1)-1, nx-p_rect_sel(2, 1), p_rect_sel(1, 2)-1, ny-p_rect_sel(2, 2)];
        bd = max(bd, bd_user);
    end

    function cb_crt_txy(opt)
        sigma_gt = str2double(ui_edt_sigma.String);
        
        if(sigma_gt>sigma_g_max)
            return
        end

        fcn_msn_start();
        
        idx_ref = ui_pu_idx_ref.Value;
        at_p(:, 5:6) = at_p(:, 5:6) - at_p(idx_ref, 5:6);
        
        rg_n_it = str2num(ui_edt_rg_n_it.String); %#ok<ST2NM>
        b_fit = true;
        sigma_rt = ilm_sigma_g_2_sigma_r(nx, ny, sigma_gt);
        radius = sigma_rt;
        p = 0.05;
        
        % calculate shifting between images
        st_opt = fcn_st_preproc_st_opt();
        bd = fcn_atp_p_rect_sel_2_bd(at_p, p_rect_sel);
        at_p_r = ilm_fd_tr_2d_bi(system_config, data, p, sigma_gt, radius, bd, at_p, rg_n_it, b_fit, st_opt);
        
		% center shifts
        if idx_ref>0
            if(opt==1)
                at_p(:, 5) = at_p_r(:, 5) - at_p_r(idx_ref, 5);
            elseif(opt==2)
                at_p(:, 6) = at_p_r(:, 6) - at_p_r(idx_ref, 6);
            else
                at_p(:, 5:6) = at_p_r(:, 5:6) - at_p_r(idx_ref, 5:6);
            end
        else
            at_p = ilm_center_tr_2d_using_at_v(at_p);
        end

        % fill in ui_lb_data
        fcn_set_at_p_2_ui_lb_data(at_p);
        
        % plot fimage
        fcn_show_ax_f_1(ui_cb_rcfft.Value);
        
        % select item in ui_lb_data
        cb_show_sel_ui_lb_data();
        
        fcn_msn_end();
    end

    function cb_crt_A_txy(source, event)
        sigma_gt = str2double(ui_edt_sigma.String);
        
        if(sigma_gt>sigma_g_max)
            return
        end

        fcn_msn_start();
        
        rg_n_it = str2num(ui_edt_rg_n_it.String); %#ok<ST2NM>
        
        % find affine transformation
        at_p = ilm_run_rg_stem(system_config, data, p, sigma_gt, at_p, rg_n_it);

        % fill in ui_lb_data
        fcn_set_at_p_2_ui_lb_data(at_p);
        
        % plot fimage
        fcn_show_ax_f_1(ui_cb_rcfft.Value);
        
        % show pcf
        %fcn_show_rg_data(sigma_gt);

        % select item in ui_lb_data
        cb_show_sel_ui_lb_data();
        
        fcn_msn_end();
    end

    function cb_crt_nr_grid(source, event)
        fcn_msn_start();
        
        sigma_gt = str2double(ui_edt_sigma.String);
        sigma_r = 1/(2*pi*sigma_gt*min(1/nx, 1/ny));
        alpha = str2double(ui_edt_alpha.String);
        
        nrg_n_it_i = str2num(ui_edt_nrg_n_it_i.String); %#ok<ST2NM>
        nrg_n_it_o = str2num(ui_edt_nrg_n_it_o.String); %#ok<ST2NM>
        
        opt = ui_pu_opt.Value;
        
        bd = fcn_atp_p_rect_sel_2_bd(at_p, p_rect_sel);
        nr_grid = ilm_run_nr_stem(system_config, data, bd, alpha, sigma_r, nr_grid, nrg_n_it_i, nrg_n_it_o, opt);
        
        fcn_msn_end();
    end

    function fcn_draw_circle_ax_f_1(bb_mouse_click)
        C = get(ax_f(1), 'CurrentPoint');
        sigma_gt = sqrt((C(1, 1)-fs_x_c)^2+(C(1, 2)-fs_y_c)^2);
        sigma_rt = ilm_sigma_g_2_sigma_r(nx, ny, sigma_gt);

        children = get(ax_f(1), 'Children');
        if bb_mouse_click
            delete(findobj(children, 'Type', 'Line'));
            
             if(sigma_gt>sigma_g_max)
                return
            end
        else
            delete(findobj(children, 'Type', 'Line', '-and', 'Tag', 'm'));
        end

        axes(ax_f(1));
        title(ax_f(1), ['sigma_g = ', num2str(sigma_gt, '%6.2f'), ' px. , sigma_r = ', num2str(sigma_rt, '%6.2f'), ' px.']);

        if bb_mouse_click
            ui_edt_sigma.String = num2str(sigma_gt, '%5.2f');
                
            ilm_plot_circle(fs_x_c, fs_y_c, sigma_gt, 'red', 'f');
            
            cb_show_sel_ui_lb_data();
        else
            ilm_plot_circle(fs_x_c, fs_y_c, sigma_gt, 'green', 'm');
        end
    end
 
    function fcn_draw_rectangle_ax_f_2(bb_mouse_click)
        C = get(ax_f(2), 'CurrentPoint');
        C = [max(1, min(nx, round(C(1, 1)))), max(1, min(ny, round(C(1, 2))))];

        d = abs(p_rect_sel-C);
        [~, ii] = min(d(:));

        if(bb_ctrl_kpf)
            if(ii<=2)
                bx = ilm_ifelse(ii==1, C(1), nx-C(1));
                bx = round((nx-ilm_pn_fact(nx-2*bx, 4))/2);
                ix_0 = bx+1;
                ix_e = nx-bx;
                p_rect_sel(:, 1) = [ix_0;ix_e];
            else
                by = ilm_ifelse(ii==3, C(2), ny-C(2));
                by = round((ny-ilm_pn_fact(ny-2*by, 4))/2);
                iy_0 = by+1;
                iy_e = ny-by;
                p_rect_sel(:, 2) = [iy_0;iy_e];
            end
        elseif(bb_shift_kpf)
            if(ii==1)
                ix_e = p_rect_sel(2, 1);
                ix_0 = round(ix_e - ilm_pn_fact(round(ix_e-C(1)+1))+1);
                p_rect_sel(1, 1) = ilm_set_bound(ix_0, 1, nx);
            elseif(ii==2)
                ix_0 = p_rect_sel(1, 1);
                ix_e = round(ilm_pn_fact(round(C(1)-ix_0+1))+ ix_0-1);
                p_rect_sel(2, 1) = ilm_set_bound(ix_e, 1, nx);
            elseif(ii==3)
                iy_e = p_rect_sel(2, 2);
                iy_0 = round(iy_e - ilm_pn_fact(round(iy_e-C(2)+1))+1);
                p_rect_sel(1, 2) = ilm_set_bound(iy_0, 1, ny);
            else
                iy_0 = p_rect_sel(1, 2);
                iy_e = round(ilm_pn_fact(round(C(2)-iy_0+1))+ iy_0-1);
                p_rect_sel(2, 2) = ilm_set_bound(iy_e, 1, ny);
            end
        end
        
        children = get(ax_f(2), 'Children');
        if bb_mouse_click
            delete(findobj(children, 'Type', 'Line'));
        else
            delete(findobj(children, 'Type', 'Line', '-and', 'Tag', 'm'));
        end
        
        axes(ax_f(2));
        if bb_mouse_click
            cb_show_sel_ui_lb_data();
        else
            ilm_plot_rectangle(p_rect_sel, 'green', 'm');
        end

        ui_pb_data_info.String = fcn_str_data_info(p_rect_sel);
    end     

    function cb_mouse_move(source, event) 
        if (bb_shift_kpf || bb_ctrl_kpf)
            b_ax_f = ilm_gui_is_over_object(source, ax_f);

            if(b_ax_f(1))
                if(bb_shift_kpf)
                    fcn_draw_circle_ax_f_1(false);
                    children = get(ax_f(2), 'Children');
                    delete(findobj(children, 'Type', 'Line', '-and', 'Tag', 'm'));
                end
            elseif(b_ax_f(2))
                if((bb_shift_kpf || bb_ctrl_kpf)&& ui_rb_select.Value)
                    fcn_draw_rectangle_ax_f_2(false);
                    children = get(ax_f(1), 'Children');
                    delete(findobj(children, 'Type', 'Line', '-and', 'Tag', 'm'));
                end
            end
        end
    end

    function cb_mouse_click(source, event)
        if (bb_shift_kpf || bb_ctrl_kpf)
            b_ax_f = ilm_gui_is_over_object(source, ax_f);

            if(b_ax_f(1))
                if(bb_shift_kpf)
                    fcn_draw_circle_ax_f_1(true);
                end
            elseif((bb_shift_kpf || bb_ctrl_kpf)&& ui_rb_select.Value)
                fcn_draw_rectangle_ax_f_2(true);
            end
        end
    end

    function cb_press_key(obj, event)
        if(strcmpi(event.Key, 'delete'))
            b_ui_lb_data = ilm_gui_is_over_object(obj, ui_lb_data);
            if(~b_ui_lb_data)
                return
            end  

            fcn_msn_start();
            
            ik_del = ui_lb_data.UserData{ui_lb_data.Value};
            data(:, :, ik_del) = [];
            n_data = ilm_nz_data(data);
            
            at_p(ik_del, :) = [];
            nr_grid.x(:, :, ik_del) = [];
            nr_grid.y(:, :, ik_del) = [];

            % fill in ui_lb_data
            fcn_set_at_p_2_ui_lb_data(at_p);
            ui_lb_data.Value = min(max(1, ik_del), n_data);
            
            % select item in ui_lb_data
            cb_show_sel_ui_lb_data();
            
            % plot fimage
            fcn_show_ax_f_1(ui_cb_rcfft.Value);
            
            fcn_msn_end();
        end
        
        bb_shift_kpf = strcmpi(event.Modifier, 'shift');
        bb_shift_kpf = ilm_ifelse(isempty(bb_shift_kpf), 0, bb_shift_kpf);
        
        bb_ctrl_kpf = strcmpi(event.Modifier, 'control');
        bb_ctrl_kpf = ilm_ifelse(isempty(bb_ctrl_kpf), 0, bb_ctrl_kpf);
    end

    function cb_release_key(source, event)
        bb_shift_kpf = 0;
        bb_ctrl_kpf = 0;
    end

    function cb_zoom(source, event)       
        fz = 1;
        if(event.VerticalScrollCount<0)
            fz = 1.5;
        elseif(event.VerticalScrollCount>0)
            fz = 1/1.5;
        end
        
        b_ax_f = ilm_gui_is_over_object(source, ax_f);

        if(b_ax_f(1))
            fz_t(1) = fz_t(1)*fz;
            axes(ax_f(1));
            zoom(fz);
        elseif(b_ax_f(2))   
            fz_t(2) = fz_t(2)*fz;
            axes(ax_f(2));
            zoom(fz);
        elseif(b_ax_f(3))  
            fz_t(3) = fz_t(3)*fz;
            axes(ax_f(3));
            zoom(fz);
        end
    end 

    function cb_close(source, event)
        b_stop = true;

        if(~bb_nr)
            nr_grid = ilm_at_p_2_nr_grid(at_p, nx, ny);
        end
        
        delete(fh);
    end
end