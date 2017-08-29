function [] = MEX(option, m_file, src, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Format input data %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(strcmp('../', src(1:3)))
        [pathstr,~,~] = fileparts(pwd);
        src = [pathstr, filesep, src(4:end)];
    end
    
    nVarargs = length(varargin);
    for k = 1:nVarargs
        varargin{k} = strcat(src, filesep, char(varargin{k}));
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% set cuda path %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    if(~exist(CUDA_PATH,'dir'))
        disp(['Cuda path not found'])
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
    
    %%%%%%%%%%%%%%% set NVCC compiler location %%%%%%%%%%%%%%%%%%
    setenv('CUDA_PATH', CUDA_PATH);    
    setenv('CUDA_BIN_PATH', CUDA_BIN_PATH);
    setenv('CUDA_LIB_PATH', CUDA_LIB_PATH);
    setenv('MW_NVCC_PATH', CUDA_BIN_PATH);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% https://developer.nvidia.com/cuda-gpus %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       %%%%%%%%%%%%%%%% set card architecture %%%%%%%%%%%%%%%%%%
    SINGLE_CARD='-gencode=arch=compute_60,code=&#92;&quot;sm_60,compute_60&#92;&quot;';
    MULTI_CARD=['-gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30'...
        ' -gencode=arch=compute_35,code=&#92;&quot;sm_35,compute_35&#92;&quot;'...
        ' -gencode=arch=compute_50,code=&#92;&quot;sm_50,compute_50&#92;&quot;'...
        ' -gencode=arch=compute_60,code=&#92;&quot;sm_60,compute_60&#92;&quot;'];		

    ARCH_FLAGS = MULTI_CARD;

    if(ispc)
        mex_cuda_filename = 'mex_CUDA_win64.xml';
    elseif(ismac)
        mex_cuda_filename = 'mex_CUDA_maci64.xml';
    else
        mex_cuda_filename = 'mex_CUDA_glnxa64.xml';  
    end 
    
       %%%%%%%%%%%%% read template mex_cuda file %%%%%%%%%%%%%%%
    mex_cuda_file_in = [pathstr, filesep, 'matlab_functions', filesep, mex_cuda_filename];
    mex_cuda = fileread(mex_cuda_file_in);
    mex_cuda = strrep(mex_cuda, 'XXX_ARCH_FLAGS', ARCH_FLAGS);
        
       %%%%%%%%%%%%%%%%% save mex_cuda file %%%%%%%%%%%%%%%%%%%%
    mex_cuda_file_out = [pwd, filesep, mex_cuda_filename];
    fid = fopen(mex_cuda_file_out, 'wt');
    fwrite(fid, mex_cuda);
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% set library path %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if(ispc)
        % OS = Windows 10
        FFTW_LIB_PATH = src;
        FFTW_LIBS = '-lfftw3f-3 -lfftw3-3';
        BLAS_LAPACK_LIB_PATH = src;
        BLAS_LAPACK_LIBS = '-lblas -llapack';     
    elseif(ismac)
        % OS = Windows 10
        FFTW_LIB_PATH = [];
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        BLAS_LAPACK_LIB_PATH = [];
        BLAS_LAPACK_LIBS = '-lblas -llapack';
    else
        % OS = scientific linux
%         FFTW_LIB = '-/opt/local/lib';  
%         FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
%         LAPACK_LIB = '/opt/local/lib';
%         BLAS_LAPACK_LIBS = '-lblas -llapack';

        % OS = ubuntu
        FFTW_LIB_PATH = '/usr/lib/x86_64-linux-gnu';  
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        BLAS_LAPACK_LIB_PATH = '/usr/lib';
        BLAS_LAPACK_LIBS = '-lblas -llapack';
    end
    
    if(~exist(FFTW_LIB_PATH,'dir'))
        disp(['FFTW path not found'])
    end
    
    if(~exist(BLAS_LAPACK_LIB_PATH,'dir'))
        disp(['Blas and Lapack paths not found'])
    end
    
    INCLUDE = ['-I"',CUDA_INC_PATH,'" -I"', src, '"'];
    LIBRARY = ['-L"',CUDA_LIB_PATH,'" ',CUDA_LIBS,' -L"',FFTW_LIB_PATH,'" ',FFTW_LIBS,' -L"',BLAS_LAPACK_LIB_PATH,'" ',BLAS_LAPACK_LIBS];

    if (strcmpi(option, 'release'))
        mex_comand = 'mex -v -largeArrayDims';
    else
        mex_comand = 'mex -v -g -largeArrayDims'; 
    end

    OUTDIR = ['-outdir ..', filesep, 'mex_bin'];

    textcommands = strjoin({mex_comand, OUTDIR, INCLUDE, LIBRARY, m_file, strjoin(varargin)});
    disp(textcommands);
    eval(textcommands);
    
    delete(mex_cuda_file_out);
end