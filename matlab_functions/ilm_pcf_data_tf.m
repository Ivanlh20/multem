function [im_rp, im_sp] = ilm_pcf_data_tf(im_r, im_s, opt, sigma_f)
    if nargin<3
        opt=1;
    end
    
    if nargin<4
        sigma_f = 4;
    end
    
    im_rp = double(im_r);
    im_sp = double(im_s);
    
    if opt==2
        im_rp = del2(im_rp);
        im_sp = del2(im_sp);
    elseif opt==3
        im_rp = ilm_lcwt(im_rp, sigma_f);
        im_sp = ilm_lcwt(im_sp, sigma_f); 
    else
        im_rp = im_rp - mean(im_rp, 'all');
        im_sp = im_sp - mean(im_sp, 'all');  
    end
    
end