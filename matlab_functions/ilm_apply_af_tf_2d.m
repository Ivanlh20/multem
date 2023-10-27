% apply affine transformation
function[im_o] = ilm_apply_af_tf_2d(system_config, im_i, dx, dy, A, txy, bg)
    [ny, nx] = size(im_i, [1, 2]);
    bsx = nx*dx;
    bsy = ny*dy;
    
    A = inv(A);
    txy = -1*A*reshape(txy, [], 1);
    
    [rx, ry, ~] = ilm_rs_grid_2d(nx, ny, bsx, bsy, 0);
    p = [rx(:), ry(:)];
    p = p*(A.') + txy.';
    rx_o = reshape(p(:, 1), [ny, nx]);
    ry_o = reshape(p(:, 2), [ny, nx]);
    
    
%     im_o = griddata(rx, ry, im_i, rx_o, ry_o, 'natural');
    im_o = interp2(rx, ry, im_i, rx_o(:), ry_o(:), 'cubic', bg);
    im_o = reshape(im_o, [ny, nx]);
    
    im_o = max(min(im_i(:)), im_o);
    
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