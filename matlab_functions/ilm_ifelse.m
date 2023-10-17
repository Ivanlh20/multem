function [c] = ilm_ifelse(bb, a, b)
	if(bb) %#ok<ALIGN>
        c = a;
    else
        c = b;
    end
end