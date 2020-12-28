function [] = ilm_mex(option, m_file, src, varargin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Format input data %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(strcmp('../', src(1:3)))
        [pathstr, ~, ~] = fileparts(pwd);
        src = [pathstr, filesep, src(4:end)];
    else
        [pathstr, ~, ~] = fileparts(src);
    end
    
    nVarargs = length(varargin);
    for k = 1:nVarargs
        varargin{k} = strcat(src, filesep, char(varargin{k}));
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% set cuda path %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CUDA_PATH = getenv('CUDA_PATH');
    
    CUDA_PATH = '';
    
    if(isempty(CUDA_PATH))
        if(ispc)
            CUDA_PATH = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.0';
        elseif(ismac)
            CUDA_PATH = '/Developer/NVIDIA/CUDA-10.0';
        else
            CUDA_PATH = '/usr/local/cuda-10.0';
        end
    end

    CUDA_PATH = ilm_replace_filesep(CUDA_PATH);
    
    %%%%%%%%%%%%%%%%%%%% get cuda version  %%%%%%%%%%%%%%%%%%%%%%%
    idx = strfind(CUDA_PATH, filesep);
    CUDA_VERSION = CUDA_PATH((idx(end)+1):end);
    CUDA_VERSION = CUDA_VERSION(~isletter(CUDA_VERSION));
    CUDA_VERSION = erase(CUDA_VERSION, {'-', '_'});
    
    if isempty(CUDA_VERSION)
        CUDA_VERSION_PATH = [CUDA_PATH, filesep,'version.txt'];
        FID = fopen(CUDA_VERSION_PATH);
        data = textscan(FID,'%s');
        fclose(FID);
        stringData = string(data{:});
        stringData = strsplit(stringData{3}, '.');
        CUDA_VERSION = [stringData{1}, '.', stringData{2}];
    end
    
    if(~exist(CUDA_PATH, 'dir'))
        disp('Error - Cuda path not found')
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
    setenv('CUDA_NVVM_PATH', [CUDA_PATH, filesep, 'nvvm']);
    setenv('MW_NVCC_PATH', CUDA_BIN_PATH);
    
    %%%%%%%%%%%%%%%%%%% set card architecture %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% https://developer.nvidia.com/cuda-gpus %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CARD_30="-gencode=arch=compute_30,code=sm_30";
    CARD_35="-gencode=arch=compute_35,code=&#92;&quot;sm_35,compute_35&#92;&quot;";
    CARD_50="-gencode=arch=compute_50,code=&#92;&quot;sm_50,compute_50&#92;&quot;";
    CARD_60="-gencode=arch=compute_60,code=&#92;&quot;sm_60,compute_60&#92;&quot;";
    CARD_70="-gencode=arch=compute_70,code=&#92;&quot;sm_70,compute_70&#92;&quot;";
    CARD_MULT = join([CARD_30, CARD_35, CARD_50, CARD_60, CARD_70], ' ');

    if 1
        gpu = gpuDevice();
        gpu_card = round(str2double(gpu.ComputeCapability));

        switch gpu_card
            case 3
                ARCH_FLAGS = CARD_30;
            case 5
                ARCH_FLAGS = CARD_50;
            case 6
                ARCH_FLAGS = CARD_60;
            case 7
                ARCH_FLAGS = CARD_70;
            otherwise
                ARCH_FLAGS = CARD_30;
        end
    else
        ARCH_FLAGS = CARD_MULT;
    end
    
    %%%%%%%%%%%%%%% read template mex_cuda file %%%%%%%%%%%%%%%%%
    if(ispc)
        mex_cuda_filename = 'mex_CUDA_win64.xml';
    elseif(ismac)
        mex_cuda_filename = 'mex_CUDA_maci64.xml';
    else
        mex_cuda_filename = 'mex_CUDA_glnxa64.xml';
    end 
    
    %%%%%%%%%%%%%%%%%% read mex_cuda template %%%%%%%%%%%%%%%%%%%%
    mex_cuda_file_in = [pathstr, filesep, 'matlab_functions', filesep, mex_cuda_filename];
    mex_cuda = fileread(mex_cuda_file_in);
    
    %%%%%%%%%%%%% replace string in mex_cuda file %%%%%%%%%%%%%%%%
    mex_cuda = strrep(mex_cuda, 'XXX_CUDA_VER', CUDA_VERSION);
    mex_cuda = strrep(mex_cuda, 'XXX_ARCH_FLAGS', ARCH_FLAGS);
    
    %%%%%%%%%%%%%%%%%%% save mex_cuda file %%%%%%%%%%%%%%%%%%%%%%
    mex_cuda_file_out = [pwd, filesep, mex_cuda_filename];
    fid = fopen(mex_cuda_file_out, 'w+');
    fwrite(fid, mex_cuda);
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% set library path %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(ispc)
        % OS = Windows 10
        FFTW_LIB_PATH = src;
        FFTW_LIBS = '-lfftw3f-3 -lfftw3-3';
        
        BLAS_LIB_PATH = src;
        BLAS_LIBS = '-lblas';
        
        LAPACK_LIB_PATH = src;
        LAPACK_LIBS = '-llapack';
        
    elseif(ismac)
        % OS = /usr/bin/ld: 
        FFTW_LIB_PATH = [];
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        BLAS_LIB_PATH = [];
        BLAS_LIBS = '-lblas';
        
        LAPACK_LIB_PATH = [];
        LAPACK_LIBS = '-llapack';
    else
        % OS = scientific linux
        %FFTW_LIB = '-/opt/local/lib';
        %FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
        %LAPACK_LIB = '/opt/local/lib';
        %BLAS_LAPACK_LIBS = '-lblas -llapack';

        % OS = ubuntu
        FFTW_LIB_PATH = '/usr/lib/x86_64-linux-gnu';
        FFTW_LIBS = '-lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
		
        BLAS_LIB_PATH = '/usr/lib/x86_64-linux-gnu/blas';
        BLAS_LIBS = '-lblas';
		
        LAPACK_LIB_PATH = '/usr/lib/x86_64-linux-gnu/lapack';
        LAPACK_LIBS = '-llapack';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BLAS_LIBS = '-lmwblas';
    LAPACK_LIBS = '-lmwlapack ';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(~exist(FFTW_LIB_PATH, 'dir'))
        disp('FFTW path not found')
    end
    
    if(~exist(BLAS_LIB_PATH, 'dir'))
        disp('Blas path not found')
    end
    
    if(~exist(LAPACK_LIB_PATH, 'dir'))
        disp('Lapack path not found')
    end
    
%     if(~ispc && ~ismac)
%         mex_cuda_filename = 'nvcc_g++.xml';
%     end 
    
    INCLUDE = ['-I"', CUDA_INC_PATH, '" -I"', src, '"'];
    LIBRARY = ['-L"', CUDA_LIB_PATH, '" ', CUDA_LIBS, ' -L"', FFTW_LIB_PATH, '" ', FFTW_LIBS, ' -L"', BLAS_LIB_PATH, '" ', BLAS_LIBS, ' -L"', LAPACK_LIB_PATH, '" ', LAPACK_LIBS];    
    mex_comand = ['mex -R2017b -f ', mex_cuda_filename, ' -v'];
    
    if (strcmpi(option, 'debug'))
        mex_comand = [mex_comand, ' -g'];
    end

    OUTDIR = ['-outdir ..', filesep, 'mex_bin'];

    textcommands = strjoin({mex_comand, OUTDIR, INCLUDE, LIBRARY, m_file, strjoin(varargin)});
    disp(textcommands);
    eval(textcommands);

    delete(mex_cuda_file_out);
end