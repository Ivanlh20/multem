function [] = MEX(option, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
	varargin{k} = strcat(path, filesep, char(varargin{k}));
end;

CUDA_PATH = getenv('CUDA_PATH');
if(isempty(CUDA_PATH))
	CUDA_INC='-I/opt/cuda/include';
else
	CUDA_INC = strcat('-I"', CUDA_PATH, filesep, 'include"');
end;

if(ispc)
    FFTW_LIB = ['-L' path];
    FFTW_LIB = [FFTW_LIB ' -lfftw3f -lfftw3'];
elseif(ismac)
    FFTW_LIB = ' -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads';
else
    FFTW_LIB = '-L/opt/local/lib'; %scientific linux
    FFTW_LIB = [FFTW_LIB ' -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads'];    
end;

ADD_INC = strcat('-I', path);
OUTDIR = strcat('..', filesep, 'mex_executables');

if (strcmpi(option, 'debug'))
    mex_comand = 'mex -v -g -largeArrayDims -outdir';
else
    mex_comand = 'mex -v -largeArrayDims -outdir';   
end;

textcommands = strjoin({mex_comand, OUTDIR, ADD_INC, CUDA_INC, mfile, strjoin(varargin), FFTW_LIB});  
disp(textcommands);
eval(textcommands);