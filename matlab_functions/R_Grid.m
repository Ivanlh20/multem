function [Rx, Ry, R2] = R_Grid(nx, ny, lx, ly, sft)
if(nargin==4)
    sft = 1;
end
if(nargin==2)
    lx = 1;
    ly = 1;
    sft = 1;
end
dRx = lx/nx; dRy = ly/ny;
Rxl = (0:1:(nx-1))*dRx;
Ryl = (0:1:(ny-1))*dRy;
[Rx, Ry] = meshgrid(Rxl, Ryl);
if(sft==1)
    Rx = ifftshift(Rx);
    Ry = ifftshift(Ry);
end
R2 = Rx.^2 + Ry.^2;