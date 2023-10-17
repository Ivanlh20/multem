function [pxy] = ilm_fit_pcf_2d(pcf_2d, dx, dy, sigma_g_px)
    [ny, nx] = size(pcf_2d, [1, 2]);
    bsx = nx*dx;
    bsy = ny*dy;
    sigma_g = sigma_g_px*min(1/bsx, 1/bsy);
    
    [~, ixy] = max(pcf_2d(:));
    [I_row, I_col] = ind2sub([ny, nx], ixy);
    pxy = [I_col, I_row].*[dx, dy];
    sigma_r = 1/(2*pi*sigma_g);
    radius = 1.5*sigma_r;
    
    pxy = ilc_fit_peak_pos_2d(pcf_2d, dx, dy, pxy, sigma_r, radius);
end