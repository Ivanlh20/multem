function [] = MEX(option, m_file, src, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
	varargin{k} = strcat(src, filesep, char(varargin{k}));
end

CUDA_PATH = getenv('CUDA_PATH');
if(isempty(CUDA_PATH))
    if(ispc)
        setenv('CUDA_PATH', 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0');
    elseif(ismac)
        setenv('CUDA_PATH', '/Developer/NVIDIA/CUDA-8.0');
    else
        setenv('CUDA_PATH', '/usr/local/cuda-8.0');    
    end    
end
CUDA_INC_BIN = [CUDA_PATH, filesep, 'bin'];
CUDA_INC_PATH = [CUDA_PATH, filesep, 'include'];

% NVCC compiler location
setenv('MW_NVCC_PATH', CUDA_INC_BIN);

if(ispc)
    FFTW_LIB = ['-L"',src,'" -lfftw3f-3 -lfftw3-3'];
    LAPACK_LIB = ['-lblas -llapack'];
    CUDA_LIB_PATH  = [CUDA_PATH, filesep, 'lib', filesep, 'x64'];      
elseif(ismac)
    FFTW_LIB = ' -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
    LAPACK_LIB = ' -lblas -llapack';
    CUDA_LIB_PATH  = [CUDA_PATH, filesep, 'lib64'];
else
    % scientific linux
%     FFTW_LIB = '-L"/opt/local/lib" -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';  
%     LAPACK_LIB = '-L"/opt/local/lib" -lblas -llapack';
    % ubuntu
    FFTW_LIB = '-L"/usr/lib/x86_64-linux-gnu" -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';  
    LAPACK_LIB = '-L"/usr/lib" -lblas -llapack';
    CUDA_LIB_PATH  = [CUDA_PATH, filesep, 'lib64'];
end

INCLUDE = ['-I"',CUDA_INC_PATH,'" -I"', src, '"'];
LIBRARY = ['-L"',CUDA_LIB_PATH,'" -lcudart -lcufft -lcublas ', FFTW_LIB, ' ', LAPACK_LIB];

if (strcmpi(option, 'release'))
    mex_comand = 'mex -v -largeArrayDims';
else
    mex_comand = 'mex -v -g -largeArrayDims'; 
end

OUTDIR = ['-outdir ..', filesep, 'mex_bin'];

textcommands = strjoin({mex_comand, OUTDIR, INCLUDE, LIBRARY, m_file, strjoin(varargin)});
disp(textcommands);
eval(textcommands);