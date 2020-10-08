classdef detector
% Detector object for Multem Input
    properties
        % STEM Detector 
        % detector_type: eDT_Circular = 1, eDT_Radial = 2, eDT_Matrix = 3
        type(1,1) uint64 {mustBeLessThanOrEqual(type,3),mustBePositive} = 1;
        cir = struct('inner_ang', 60, 'outer_ang', 180);                        % Inner & Outer angle(mrad) 
        radial = struct('x', 0, 'fx', 0);                                       % radial detector angle(mrad) & radial sensitivity value
        matrix = struct('R', 0, 'fR', 0);                                       % 2D detector angle(mrad) & 2D sensitivity value
    end
end