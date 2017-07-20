function [x] = ilm_rand(x_min, x_max)
    x = x_min+(x_max-x_min)*rand();
end