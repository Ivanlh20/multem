function[system_config] = ilm_dflt_system_conf(device, id_gpu)
    if(nargin<2)
        id_gpu = 0;
    end
	
    if(nargin<1)
        device = 2;
    end   
	
    system_config.device = device; % eD_CPU = 1, eD_GPU = 2
    system_config.precision = 1; % eP_Float = 1, eP_double = 2
    
    %%%%%%%%%%%%%%%%%%%%%%%% cpu config %%%%%%%%%%%%%%%%%%%%%%%%
    system_config.cpu_n_proc = 1;
    system_config.cpu_n_thread = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%% gpu config %%%%%%%%%%%%%%%%%%%%%%%%
    system_config.gpu_device = id_gpu;
%  system_config.gpu_n_stream = 1;
%  system_config.gpu_device = 0;
end
