function [im_o] = ilm_minmax_norm_image(im_i)
    [im_min, im_max] = ilm_min_max(im_i);
    im_o = (im_i-im_min)/(im_max-im_min);
end