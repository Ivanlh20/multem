function xo = tfm_pn_fact(x, b_min)
    % Returns an integer corresponding to the closest prime number
    % factorisation for computational efficiency (mainly relevant for CUDA)
    % Inputs:
    %   x       Input number
    %   b_min   (optional) either: 
    %               1 (nearest solution) (default)
    %               2 (next smaller solution)
    %               3 (next larger solution)
    
    x = int64(x);
    x = min(max(2,x),2^16);
    p_min = 64;
    
    if nargin < 2
        b_min = 1;
    end
    
    switch b_min
        case 1
            x_t(1) = fcn_pn_down(x);
            x_t(2) = fcn_pn_up(x);
            [~,id] = sort(abs(x-x_t));
            xo = x_t(id(1));
        case 2
            x = max(3,x);
            xo = fcn_pn_down(x-1);
        case 3
            xo = fcn_pn_up(x+1);
    end
    xo = double(xo);
    function x = fcn_pn_up(x)
        pn = factor(x);
        b_pn = fcn_b_pn(x,pn);
        while b_pn
            x = x + 1;
            pn = factor(x);
            b_pn = fcn_b_pn(x,pn);
        end
    end
    
    function x = fcn_pn_down(x)
        pn = factor(x);
        b_pn = fcn_b_pn(x,pn);
        while b_pn
            x = x - 1;
            pn = factor(x);
            b_pn = fcn_b_pn(x,pn);
        end
    end

    function b_pn = fcn_b_pn(x,pn)
        if x > p_min
                b_pn = any(pn>7) || ~any(pn==2);
            else
                b_pn = ~all(pn==2);
        end
    end
end