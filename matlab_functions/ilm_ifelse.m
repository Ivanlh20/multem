function [c] = ilm_ifelse(bb, a, b)
	if(bb)
        c = a;
    else
        c = b;
    end
end