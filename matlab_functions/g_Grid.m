function [gx, gy, g2, gmax] = g_Grid(nx, ny, lx, ly, sft)
if(nargin==4)
    sft = 1;
end
if(nargin==2)
    lx = 1;
    ly = 1;
    sft = 1;
end
nxh = nx/2; nyh = ny/2;
dgx = 1/lx; dgy = 1/ly;
gxl = (-nxh:1:(nxh-1))*dgx;
gyl = (-nyh:1:(nyh-1))*dgy;
[gx, gy] = meshgrid(gxl, gyl);
if(sft==1)
    gx = ifftshift(gx);
    gy = ifftshift(gy);
end
g2 = gx.^2 + gy.^2;
gmax = (2/3)*min([nxh*dgx, nyh*dgy]);