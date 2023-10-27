function[mtf] = ilm_create_psf_g(coef_mtf, nx, ny, shift_bb, support_bb)
%     It need to be finished
    if nargin<4
        support_bb = true;
    end

    if nargin<4
        shift_bb = false;
    end

    nxh = nx/2;
    nyh = ny/2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The parameterization is given in fractions of Nyquist frequency = 1/(2*dx)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % g_nyquist = 1/(2*dx) = nx/(2*lx) = nx*dgx/2 = 1 --> dgx = 2/nx;
    % 
    drx = 0.5;
    dry = 0.5;
    [rx, ry] = meshgrid((-nxh:(nxh-1))*drx, (-nyh:(nyh-1))*dry);
    r = sqrt(rx.^2+ry.^2);
    
    % matlab sinc function is defines as sin(pi x)/(pi x)
%     mtf = mtf_par.fun(g);
    mtf_rad = ilm_eval_rad_mtf(coef_mtf, r);
    
    if support_bb
        mtf = mtf_rad.*sinc(rx/2).*sinc(ry/2);
    else
        mtf = mtf_rad;
    end
    
    if shift_bb
        mtf = fftshift(mtf);
    end
end