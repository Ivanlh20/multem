function[h, h_min, h_max, dh] = ilm_hist(x, n_bin, n_left_space, n_right_space)
    if nargin<4
        n_right_space = 0;
    end
    
    if nargin<3
        n_left_space = 0;
    end
    
    n_bin_r = n_bin - n_left_space - n_right_space; 
    
    x = double(x);
    h_min = min(x);
    h_max = max(x);
    dh = (h_max - h_min)/n_bin_r; 
    idx = floor((x - h_min)/dh)+1;
    idx = n_left_space + min(n_bin_r, max(1, idx));
    
    h = zeros(n_bin, 1, 'double');
    n_idx = numel(idx);
    for ix=1:n_idx
        ik_s = idx(ix);
        h(ik_s) = h(ik_s) + 1;
    end
end