function [parm] = ilm_build_ellip_parm_interp2(coef, nx, ny, nx_r, ny_r, radius_r)
    p_c_r = [nx_r, ny_r]/2 + 1;
    
    p_c = coef(1:2);
    rad_a = coef(3);
    rad_b = coef(4);
    theta = coef(5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mr = ilm_rot_mx_2d(-theta);
    ms = [rad_a/radius_r, 0; 0 rad_b/radius_r];
    msc = mr*ms*(mr.');
    msc = msc.';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parm.radius_r = radius_r;
    [rx_r, ry_r] = meshgrid((1:nx_r)-p_c_r(1), (1:ny_r)-p_c_r(2));
    rxy_r = [rx_r(:), ry_r(:)]*msc;
    parm.rx_r = reshape(rxy_r(:, 1), size(rx_r));
    parm.ry_r = reshape(rxy_r(:, 2), size(ry_r));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [parm.rx, parm.ry] = meshgrid((1:nx)-p_c(1), (1:ny)-p_c(2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parm.fltr_rs = ms(1)*ms(4);
    parm.p_c_r = p_c_r;
    nxh_r = nx_r/2;
    nyh_r = ny_r/2;
    parm.ax_r = (p_c_r(1)-nxh_r):1:(p_c_r(1)+nxh_r-1);
    parm.ay_r = (p_c_r(2)-nyh_r):1:(p_c_r(2)+nyh_r-1);
    parm.nx_r = nx_r;
    parm.ny_r = ny_r;
end
