function[fbw] = ilm_wd_butwth(bs, dr, n, r_wd, sft, r_max, p_c)
    if length(p_c) == 1
        p_c = [p_c, p_c];
    end
    
    if length(r_max) == 1
        r_max = [r_max, r_max];
    end
    
    if length(r_wd) == 1
        r_wd = [r_wd, r_wd];
    end
    
    radius_x = r_wd(1);
    radius_y = r_wd(2);
    rx = (0:1:(bs(1)-1))*dr(1) - p_c(1);
    ry = ((0:1:(bs(2)-1))*dr(2)).' - p_c(2);
    
    if bs(1)>1
        fx = 1./(1 + (rx/radius_x).^(2*n));
        fx_sft = max(fx(1), fx(end));
        fx = (fx - fx_sft)/(0.5-fx_sft);
        fx = max(0, min(1.0,  fx));
    else
        fx = 1;
    end 
    
    if bs(2)>1
        fy = 1./(1 + (ry/radius_y).^(2*n));
        fy_sft = max(fy(1), fy(end));
        fy = (fy - fy_sft)/(0.5-fy_sft);
        fy = max(0, min(1.0,  fy));
    else
        fy = 1;
    end
    
    fbw = fx.*fy;
    
    if sft
       fbw = ifftshift(fbw); 
    end
end