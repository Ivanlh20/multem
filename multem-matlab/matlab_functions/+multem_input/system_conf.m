classdef system_conf
% System Configuration object for Multem Input
    properties
       %%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
       
        % eP_Float = 1, eP_double = 2
        precision(1,1) uint64 {mustBeLessThanOrEqual(precision,2),mustBePositive} = 1; 
        % eD_CPU = 1, eD_GPU = 2
        device(1,1) uint64 {mustBeLessThanOrEqual(device,2),mustBePositive} = 2;
        % Number of Cores CPU (It will be used in the future)
        cpu_ncores
        % # of CPU Threads
        cpu_nthread(1,1) uint64 {mustBePositive} = 1;
        % Select GPU (for Multi-GPU Setups)
        gpu_device(1,1) uint64 {mustBeNonnegative} = 0;
    end
end