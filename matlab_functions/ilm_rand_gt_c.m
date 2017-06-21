function [b] = ilm_rand_gt_c(a)
    a = 1-a;
	b = (rand()>a);
end