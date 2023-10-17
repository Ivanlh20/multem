function[pcf] = ilm_pcf_2d(system_config, im_r, im_s, dx, dy, sigma_g_px)
    [ny, nx] = size(im_r, [1, 2]);
    bsx = nx*dx;
    bsy = ny*dy;
    sigma_g = sigma_g_px*min(1/bsx, 1/bsy);

    [~, ~, g2] = ilm_fs_grid_2d(nx, ny, bsx, bsy, 1);
    M = exp(-0.5*g2/sigma_g^2);
    
    fpcf = conj(fft2(fftshift(im_r))).*fft2(fftshift(im_s));
    fpcf = M.*exp(1i*angle(fpcf));
    pcf = ifftshift(real(ifft2(fpcf)));
    
    if 0
        figure(2);clf;
        subplot(2, 2, 1);
        imagesc(im_r);
        axis image;
        colormap jet;
        subplot(2, 2, 2);
        imagesc(im_s);
        axis image;
        colormap jet;
        subplot(2, 2, 3);
        imagesc(M);
        axis image;
        colormap jet;
        subplot(2, 2, 4);
        imagesc(pcf);
        axis image;
        colormap jet;        
    end
end