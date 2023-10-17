function [im_o] = ilm_pcf_2d_wd_bd_tf(system_config, im_i, dx, dy, A, txy, p, bd_px)
    [ny, nx] = size(im_i, [1, 2]);
    
    im_o = ilm_apply_af_tf_2d(system_config, im_i, dx, dy, A, txy, 0);
    
    rect = ilm_bd_2_rect(nx, ny, bd_px);
    rect = max(1, round(rect));
    p_c = mean(rect, 1);
    
    radius_x = 0.5*(nx-bd_px(1)-bd_px(2))*(1-p);
    radius_y = 0.5*(ny-bd_px(3)-bd_px(4))*(1-p);
    rx = (0:1:(nx-1)) - p_c(1);
    ry = (0:1:(ny-1)).' - p_c(2);
    n = 8;
    
    fx = 1./(1 + (rx/radius_x).^(2*n));
    fx_sft = max(fx(1), fx(end));
    fx = (fx - fx_sft)/(0.5-fx_sft);
    fx = max(0, min(1.0,  fx));
    fy = 1./(1 + (ry/radius_y).^(2*n));
    fy_sft = max(fy(1), fy(end));
    fy = (fy - fy_sft)/(0.5-fy_sft);
    fy = max(0, min(1.0,  fy));
 
    im_o = fx.*fy.*im_o;
     
    if 0
        
        figure(1);clf;
        subplot(1, 2, 1);
        imagesc(im_i);
        axis image;
        colormap jet;
        subplot(1, 2, 2);
        imagesc(im_o);
        axis image;
        colormap bone;
    end
end

