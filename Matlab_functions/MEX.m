function [] = MEX(option, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
	varargin{k} = strcat(path, filesep, char(varargin{k}));
end;

CUDA_PATH = getenv('CUDA_PATH');

if(isempty(CUDA_PATH))
    if(ispc)
        setenv('CUDA_PATH', 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5');
    elseif(ismac)
        setenv('CUDA_PATH', '/Developer/NVIDIA/CUDA-7.5');
    else
%         % scientific linux
%         setenv('CUDA_PATH', '/usr/local/cuda-7.5');
        % ubuntu
        setenv('CUDA_PATH', '/usr/local/cuda-7.5');        
    end;    
end;

% NVCC compiler location
setenv('MW_NVCC_PATH', strcat(CUDA_PATH, filesep, 'bin'));
% Cuda libraries
CUDA_INC = strcat('-I"', CUDA_PATH, filesep, 'include"');

if(ispc)
    FFTW_LIB = ['-L' path ' -lfftw3f -lfftw3'];
    LAPACK_LIB = ['-L' path ' -lblas -llapack'];
elseif(ismac)
    FFTW_LIB = ' -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
    LAPACK_LIB = ' -lblas -llapack';
else
%     % scientific linux
%     FFTW_LIB = '-L/opt/local/lib -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';  
%     LAPACK_LIB = '-L/opt/local/lib -lblas -llapack';
    % ubuntu
    FFTW_LIB = '-L/usr/lib/x86_64-linux-gnu -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';  
    LAPACK_LIB = '-L/usr/lib/x86_64-linux-gnu -lblas -llapack';
end;

ADD_INC = strcat('-I', path);
OUTDIR = strcat('..', filesep, 'mex_executables');

if (strcmpi(option, 'debug'))
    mex_comand = 'mex -g -largeArrayDims -outdir';
else
    mex_comand = 'mex -largeArrayDims -outdir';   
end;

textcommands = strjoin({mex_comand, OUTDIR, ADD_INC, CUDA_INC, mfile, strjoin(varargin), LAPACK_LIB, FFTW_LIB});  
disp(textcommands);
eval(textcommands);