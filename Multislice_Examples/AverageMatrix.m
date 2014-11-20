function [Ma] = AverageMatrix(M, nx, ny, uave)

    function [Mo] = AverageInside(Mi)
        [nyi, nxi] = size(Mi);
        Mo =  zeros(nyi+1, nxi+1);
        Mo(1:nyi, 1:nxi) =  Mi;
        Mo(nyi+1, :) = Mo(1, :);
        Mo(:, nxi+1) = Mo(:, 1);
        
        M1 = flipud(Mo);
        M2 = fliplr(M1);
        M3 = flipud(M2);
        
        Mo = (Mo + M1 + M2 + M3)/4;
        Mo = Mo(1:nyi, 1:nxi);
    end

    [nyt, nxt] = size(M);
    nxc = nxt/nx; nyc = nyt/ny;
    Ma = zeros(nyc, nxc);
    for i = 1:nx
        for j = 1:ny
            Mt = M(((j-1)*nyc+1):(j*nyc), ((i-1)*nxc+1):(i*nxc));  
            if(uave)
                Ma = Ma + AverageInside(Mt);
            else
                Ma = Ma + Mt;
            end;
        end;
    end;
    Ma = Ma/(nx*ny);
end