clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
valueSet = 50:50:400;
n_E0 = length(valueSet);
keySet = 1:n_E0;
E0_map = containers.Map(keySet,valueSet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, 5f, 6d, 7p, (8s, 5g, 6f, 7d, 8p, and 9s)
% keySet = {'1s','2s','2p','3s','3p','4s','3d','4p','5s','4d','5p'};
valueSet = {'1s','2s','2p','3s','3p','4s','3d','4p','4d','5p'};
n_orb = length(valueSet);
keySet = 1:n_orb;
orb_map = containers.Map(keySet,valueSet);

% g
n_g = 2^13;
n_r_0 = 1024;
n_r = 128;
r_min = 0.02;
r_max = 16;
v_ftr_lim = 1e-4;
r_0 = fcn_log_grid(r_min, r_max, n_r_0);

% set initial paramters
load('data_table_edx.mat');

n_Z = length(data_table);
% read data
for iZ=1:n_Z
    for iorb=1:n_orb
        orb = orb_map(iorb);
        if data_table(iZ).orb(iorb).E0(1).fg(1)<1e-7
            continue;
        end
        
        for iE_0=1:n_E0
            E_0 = E0_map(iE_0);
            % read data
            g = data_table(iZ).orb(iorb).E0(iE_0).g;
            fg = data_table(iZ).orb(iorb).E0(iE_0).fg;
            
            disp([num2str(iZ), ' - ', num2str(E_0), ' - ', orb, ' - ', num2str(length(g))])

            % set fg parameters
            dg = (g(end)-g(1))/(n_g-1);
            g_i = (0:1:(n_g-1))*dg;

            if 0
                % fix maximum value at g = 0 in the raw data
                fg = fcn_fix_fg_raw(g, fg);

                % fix maximum value at g = 0 in the interpolated data
                [fg_i, bb_decay] = fcn_fix_fg_intrpl(g, fg, g_i);
            else
                fg_sft = fcn_y_sft(fg);
                fg_i = interp1(g, log(fg - fg_sft), g_i, 'linear');
                fg_i = exp(fg_i) + fg_sft;
            end

            % potential
            J_0 = besselj(0, 2*pi*g_i.*r_0);
            v_i = 2*pi*trapz(g_i, fg_i.*g_i.*J_0, 2);

            % vr smooth spline
            if 0
                v_i = fcn_vr_smooth(r_0, v_i);
            end

            % r_max
%             r_lim = fcn_get_x_lim_smooth(r_0, v_i, v_ftr_lim);
            r_lim = fcn_get_x_lim_smooth_2(r_0, v_i, v_ftr_lim);

            % set data
            vr_ss_tmp(iE_0).fg_i = fg_i;
            vr_ss_tmp(iE_0).v_i = v_i;
            vr_ss_tmp(iE_0).r_lim = r_lim;
        end

        % generate grid
        r_lim = max([vr_ss_tmp.r_lim]);
        r_i = fcn_log_grid(r_min, r_lim, n_r);

        % set data
        for iE_0=1:n_E0
            v_i = vr_ss_tmp(iE_0).v_i;
            vr_sft = fcn_y_sft(v_i);
            vr_log = log(v_i - vr_sft);
            v_log_i = spline(r_0, [0; vr_log; 0], r_i);
            v_i = exp(v_log_i) + vr_sft;
            v_i = max(v_i, 0);

            data_table(iZ).orb(iorb).E0(iE_0).r_max = vr_ss_tmp(iE_0).r_lim;
            data_table(iZ).orb(iorb).E0(iE_0).r = r_i;
            data_table(iZ).orb(iorb).E0(iE_0).vp = v_i;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            g = data_table(iZ).orb(iorb).E0(iE_0).g;
            fg = data_table(iZ).orb(iorb).E0(iE_0).fg;
            fg_i = vr_ss_tmp(iE_0).fg_i;

            fg_sft = fcn_y_sft(fg_i);
            fg_log = log(fg-fg_sft);
            fg_log_i = log(fg_i-fg_sft);

            if 0
                figure(1); clf;
                subplot(2, 2, 1);
                plot(g, fg, '-r', g_i, fg_i, '-b');
                title(['feg_edx - Z=', num2str(iZ),' - orb = ',  orb], 'interpreter', 'none');
                legend({'Tabulated', 'Ivan'});
                pbaspect([1.5 1 1])

                subplot(2, 2, 2);
                plot(r_i, v_i, '-+b');
                xlim([0 r_lim]);
                legend('Ivan');

                if 0
                    str_prefix = ['potential.E', num2str(E_0), '.', ilm_Z_2_str(iZ), '.Orbital_', orb];
                    r = eval([str_prefix, '.r_log']);
                    vr = eval([str_prefix, '.V_EDX']);
                    ftr = mean(v_i)/mean(vr);
                    plot(r, vr*ftr, '-r', r_0, v_i, '-b');
                    xlim([0 r_lim]);
                    title(['Vp_edx - Z=', num2str(iZ),' - orb = ',  orb], 'interpreter', 'none');
                    legend({'Zezhong', 'Ivan'});
                end

                pbaspect([1.5 1 1])
                subplot(2, 2, 3);
                plot(g, fg_log, '-r');
                hold on;
                plot(g_i, fg_log_i, '-b');
                title(['feg_log_edx - Z=', num2str(iZ),' - orb = ',  orb], 'interpreter', 'none');
                legend({'Tabulated', 'Ivan'});
                pbaspect([1.5 1 1])
                subplot(2, 2, 4);
                plot(r_i, v_log_i, '-b');
                title(['Vp_log_edx - Z=', num2str(iZ),' - orb = ',  orb], 'interpreter', 'none');
                pbaspect([1.5 1 1])
            end
        end
    end
end

% saving data
save('data_table_edx.mat', 'data_table');

% fix maximum value at g = 0 in the raw data
function [fg] = fcn_fix_fg_raw(g, fg)
    for ik=2:5
        if fg(ik-1)<=fg(ik)
            fg(ik) = min(fg(ik-1), fg(ik));
            g_t = [g(1:(ik-1)), g((ik+1):end)];
            fg_t = [0, fg(1:(ik-1)), fg((ik+1):end), 0];
            ss = spline(g_t, fg_t);
            fg_ik_s = ppval(ss, g(ik));
            fg(ik) = 0.5*(fg(ik) + fg_ik_s);
        end
    end
end

% fix maximum value at g = 0 in the interpolated data
function [fg_i, bb_decay] = fcn_fix_fg_intrpl(g, fg, g_i)
    bb_decay = all(diff(fg)<0);

    if bb_decay
       fg_sft = fcn_y_sft(fg);
       fg = log(fg-fg_sft);
    end
    
    g_lim = max(g)/2;
    coef = polyfit(g(g>g_lim), fg(g>g_lim), 1);
    m_e = coef(1);
    ss = spline(g, [0, fg, m_e]);
    fg_i = ppval(ss, g_i);

    [~, idx_max] = max(fg_i);
    if idx_max>1
        g_m = (g(1)+g(2))/2;
        g_t = [g(1), g_m, g(2:end)];
        while idx_max>1
            fg_m = ppval(ss, g_m);
            fg_a = fg(2) + 0.99*(fg_m-fg(2));
            fg_t = [0, [fg(1), fg_a, fg(2:end)], m_e];
            ss = spline(g_t, fg_t);
            fg_i = ppval(ss, g_i);
            [~, idx_max] = max(fg_i);
        end
    end
    
    if bb_decay
       fg_i = exp(fg_i) + fg_sft;
    end
end  
        
% vr smooth spline
function [vr_s] = fcn_vr_smooth(r, vr)
    vr_sft = fcn_y_sft(vr);
    vr_log = log(vr - vr_sft);
            
    w = ones(size(r));
    w(1) = 1e+3;
    p_a = [0.995, 0.999, 0.9999, 0.99999];
    ftr_a = [1e-5, 1e-4, 1e-3, 1e-2];
    for it_sp=1:4
        % smooth
        v_log_it = csaps(r, vr_log, p_a(it_sp), r, w);
        
        % r_max
        ir_lim = fcn_get_lim(r, exp(v_log_it) + vr_sft, ftr_a(it_sp));

        % constraint
        vr_log(ir_lim:end) = v_log_it(ir_lim:end);
    end
    
    vr_s = exp(vr_log) + vr_sft;
end
            
function [r] = fcn_log_grid(r_min, r_max, n_r)
    alpha = log(r_max/r_min+1)/(n_r-1);
    r = r_min*(exp(alpha*(0:1:(n_r-1)))-1);
    r = reshape(r, [], 1);
end

function [ix_lim, x_lim] = fcn_get_lim(x, y, ftr)
    i_min = 0.5*min(abs(y));
    i_max = max(y);
    y_lim = i_min + ftr*(i_max-i_min);
    ix_lim = find(y<y_lim);
    if isempty(ix_lim)
        ix_lim = length(y);
    end
    ix_lim = ix_lim(1);
    x_lim = x(ix_lim);
end

function x_lim = fcn_get_x_lim_smooth(x, y, ftr)
    i_min = 0*0.5*min(abs(y));
    i_max = max(y);
    y_lim = i_min + ftr*(i_max-i_min);
    x_i = x(1):1e-3:x(end);
    y_i = spline(x, y, x_i)-y_lim;
    ix_lim = find(y_i<y_lim);
    if isempty(ix_lim)
        ix_lim = length(y);
    end
    ix_lim = ix_lim(1);
    x_lim = x_i(ix_lim);
end

function x_lim = fcn_get_x_lim_smooth_2(x, y, ftr)
    i_min = 0*0.5*min(abs(y));
    i_max = max(y);
    y_lim = i_min + ftr*(i_max-i_min);
    x_i = x(1):1e-3:x(end);
    y_i = interp1(x, y, x_i, 'linear');
    ix_lim = find(y_i<y_lim);
    if isempty(ix_lim)
        x_lim = x(end);
    else
        ix_lim = ix_lim(1);
        x_lim_0 = x_i(ix_lim);

        pp = interp1(x, y, 'linear', 'pp');
        chi = @(x)ppval(pp, x)-y_lim;
        x_lim = fzero(chi, x_lim_0);
%         disp([x_lim, x_lim_0])
    end
end

function [y_sft]= fcn_y_sft(y)
    y_min = min(y);
    y_sft_p = 0;

    if y_min<=1e-12
        y_sft_p = max(1e-6, abs(y(end)));
    else
        y_min = 0;
    end
    y_sft = y_min - y_sft_p;
end