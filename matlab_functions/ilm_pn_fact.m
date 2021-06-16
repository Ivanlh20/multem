function[n_s] = ilm_pn_fact(n_r, opt)
    % eDST_Closest = 1, eDST_Less_Than = 2, eDST_Greater_Than = 3
    % eDST_ELess_Than = 4, eDST_EGreater_Than = 5
    if nargin < 2
        opt = 1;
    end
    
	b_2 = 2;
    b_3 = 3;
    b_5 = 5;
    b_7 = 7;
    
	e_2_m = 16;
    e_3_m = 7;
    e_5_m = 5;
    e_7_m = 4;

	pn_0 = b_2^5;
	pn_e = b_2^e_2_m;
    
    b_2_a = reshape(b_2.^(1:e_2_m), [], 1);
    b_3_a = reshape(b_3.^(0:e_3_m), 1, []);
    b_5_a = reshape(b_5.^(0:e_5_m), 1, 1, [], 1);
    b_7_a = reshape(b_7.^(0:e_7_m), 1, 1, 1, []);
    
    n_a = b_2_a.*b_3_a.*b_5_a.*b_7_a;
    n_a = n_a(:);
    n_a((pn_0>=n_a) | (n_a>pn_e)) = [];
    n_a = [2.^((1:5).'); unique(n_a)];
    
    [~, idx] = min(abs(n_a - n_r));
    n_s = n_a(idx);
        
    if (opt==2) && (n_s >= n_r)
       n_s = n_a(idx-1);
    end

    if (opt==3) && (n_s <= n_r)
       n_s = n_a(idx+1);
    end

    if (opt==4) && (n_s > n_r)
       n_s = n_a(idx-1);
    end
    
    if (opt==5) && (n_s < n_r)
       n_s = n_a(idx+1);
    end
end