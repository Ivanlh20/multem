function ilm_compile(option, m_file, src)
    if(strcmp('../', src(1:3)))
        src = [pwd, filesep, src(4:end)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% set cuda path %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CUDA_PATH = getenv('CUDA_PATH');
    if(isempty(CUDA_PATH))
        if(ispc)
            CUDA_PATH = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0';
        elseif(ismac)
            CUDA_PATH = '/Developer/NVIDIA/CUDA-8.0';
        else
            CUDA_PATH = '/usr/local/cuda-8.0';
        end    
    end
    
    CUDA_BIN_PATH = [CUDA_PATH, filesep, 'bin'];
    CUDA_INC_PATH = [CUDA_PATH, filesep, 'include'];
     if(ispc)
        CUDA_LIB_PATH = [CUDA_PATH, filesep, 'lib', filesep, 'x64'];
    elseif(ismac)
        CUDA_LIB_PATH = [CUDA_PATH, filesep, 'lib64'];
     else
        CUDA_LIB_PATH = [CUDA_PATH, filesep, 'lib64'];  
    end     
    CUDA_LIBS = '-lcudart -lcufft -lcublas';  
    
    % NVCC compiler location
    setenv('CUDA_PATH', CUDA_PATH);    
    setenv('CUDA_BIN_PATH', CUDA_BIN_PATH);
    setenv('CUDA_LIB_PATH', CUDA_LIB_PATH);
    setenv('MW_NVCC_PATH', CUDA_BIN_PATH);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% set library path %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if(ispc)
        FFTW_LIB_PATH = src;
        FFTW_LIBS = '-lfftw3f-3 -lfftw3-3';
        BLAS_LAPACK_LIB_PATH = src;
        BLAS_LAPACK_LIBS = '-lblas -llapack';     
    elseif(ismac)
        FFTW_LIB_PATH = [];
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        BLAS_LAPACK_LIB_PATH = [];
        BLAS_LAPACK_LIBS = '-lblas -llapack';
        
    else
        % scientific linux
%         FFTW_LIB = '-/opt/local/lib';  
%         FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
%         LAPACK_LIB = '/opt/local/lib';
%         BLAS_LAPACK_LIBS = '-lblas -llapack';

        % ubuntu
        FFTW_LIB_PATH = '/usr/lib/x86_64-linux-gnu';  
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        BLAS_LAPACK_LIB_PATH = '/usr/lib';
        BLAS_LAPACK_LIBS = '-lblas -llapack';
    end

    INCLUDE = ['-I"',CUDA_INC_PATH,'" -I"', src, '"'];
    LIBRARY = ['-L"',CUDA_LIB_PATH,'" ',CUDA_LIBS,' -L"',FFTW_LIB_PATH,'" ',FFTW_LIBS,' -L"',BLAS_LAPACK_LIB_PATH,'" ',BLAS_LAPACK_LIBS];

    if (strcmpi(option, 'release'))
        mex_comand = 'mex -v -largeArrayDims';
    else
        mex_comand = 'mex -v -g -largeArrayDims'; 
    end

    OUTDIR = ['-outdir ..', filesep, 'mex_bin'];
    
    CARD_20 = '-gencode arch=compute_20,code=sm_20';
    CARD_30 = '-gencode arch=compute_30,code=sm_30';
    CARD_35 = '-gencode arch=compute_35,code=sm_35';
    CARD_50 = '-gencode arch=compute_50,code=sm_50';
    CARD_60 = '-gencode arch=compute_60,code=sm_60';
    CARD_MULT = [CARD_20, ' ', CARD_30, ' ', CARD_35, ' ', CARD_50, ' ', CARD_60];
    ARCH_FLAGS = CARD_60;
  
    NVCC_FLAGS='-std=c++11 --default-stream legacy';
    COMP_FLAGS='--compiler-options=/Zp8,/GR,/W3,/EHs,/nologo,/MD';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% compilation %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    cmd1 = [CUDA_BIN_PATH, filsep, 'nvcc -c --compiler-options=/D_CRT_SECURE_NO_DEPRECATE,/D_SCL_SECURE_NO_DEPRECATE,/D_SECURE_SCL=0,/DMATLAB_MEX_FILE' ...
        ' ', INCLUDE,...
        ' -I"',matlabroot,'/extern/include"',...
        ' ', ARCH_FLAGS,...
        ' ', NVCC_FLAGS,...
        ' -std=c++11 --compiler-options=-ansi,-fexceptions,-fPIC,-fno-omit-frame-pointer,-lcufft,-lcufftw,-pthread -O -DNDEBUG' ...
        ' ' files{k} '.cu -o ' files{k} '.o'];
    
    cmd2 = ['/usr/bin/g++ -pthread -Wl,--no-undefined -Wl,--no-as-needed -shared -O' ...
        ' -Wl,--version-script,"/usr/local/MATLAB/R2016b/extern/lib/glnxa64/mexFunction.map"' ...
        ' ' strcat(files{k},'.o') ' -ldl' ...
        ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcusparse.so' ... % -lcusparse
        ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcublas_static.a' ... % -lcublas_static
        ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcusparse_static.a' ... % -lcusparse_static
        ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libculibos.a' ... % -lculibos'
        ' -L/usr/local/cuda-8.0/lib64 -Wl,-rpath-link,/usr/local/MATLAB/R2016b/bin/glnxa64' ...
        ' -L"/usr/local/MATLAB/R2016b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -lmwgpu -lcufft -lcufftw' ...
        ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcudart.so' ... % /usr/local/MATLAB/R2016a/bin/glnxa64/libcudart.so.7.5
        ' -o ' files{k} '.mexa64'];
    
    cmd3 = ['rm -f ' files{k} '.o'];


%     cmd1 = ['/usr/local/cuda-8.0/bin/nvcc -c --compiler-options=-D_GNU_SOURCE,-DMATLAB_MEX_FILE' ...
%         ' -I"/usr/local/cuda-8.0/include"' ...
%         ' -I"/usr/local/MATLAB/R2016b/extern/include"' ...
%         ' -I"/usr/local/MATLAB/R2016b/simulink/include"' ...
%         ' -I"/usr/local/MATLAB/R2016b/toolbox/distcomp/gpu/extern/include/"' ...
%         ' -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50' ...
%         ' -std=c++11 --compiler-options=-ansi,-fexceptions,-fPIC,-fno-omit-frame-pointer,-lcufft,-lcufftw,-pthread -O -DNDEBUG' ...
%         ' ' files{k} '.cu -o ' files{k} '.o'];
%     
%     cmd2 = ['/usr/bin/g++ -pthread -Wl,--no-undefined -Wl,--no-as-needed -shared -O' ...
%         ' -Wl,--version-script,"/usr/local/MATLAB/R2016b/extern/lib/glnxa64/mexFunction.map"' ...
%         ' ' strcat(files{k},'.o') ' -ldl' ...
%         ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcusparse.so' ... % -lcusparse
%         ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcublas_static.a' ... % -lcublas_static
%         ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcusparse_static.a' ... % -lcusparse_static
%         ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libculibos.a' ... % -lculibos'
%         ' -L/usr/local/cuda-8.0/lib64 -Wl,-rpath-link,/usr/local/MATLAB/R2016b/bin/glnxa64' ...
%         ' -L"/usr/local/MATLAB/R2016b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -lmwgpu -lcufft -lcufftw' ...
%         ' /usr/local/cuda-8.0/targets/x86_64-linux/lib/libcudart.so' ... % /usr/local/MATLAB/R2016a/bin/glnxa64/libcudart.so.7.5
%         ' -o ' files{k} '.mexa64'];
%     
%     cmd3 = ['rm -f ' files{k} '.o'];
    
    disp([files{k} '.cu'])
    if system(cmd1); error('%s failed step 1',files{k}); end
    if system(cmd2); error('%s failed step 2',files{k}); end
    if system(cmd3); error('%s failed step 3',files{k}); end
    disp('MEX completed successfully.')
end