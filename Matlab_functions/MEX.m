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

if(isunix)
    FFTW_THREAD_LIB='-lfftw3f_threads -lfftw3_threads';
else
    FFTW_THREAD_LIB='';
end;

ADD_INC = strcat('-I', path);
OUTDIR = strcat('..', filesep, 'mex_executables');

if (strcmpi(option, 'debug'))
    textcommands = strjoin({'mex -silent -g -largeArrayDims -outdir', OUTDIR, ADD_INC, CUDA_INC...
    , mfile, strjoin(varargin), ['-L' path ' -lfftw3f'], ['-L' path ' -lfftw3'], FFTW_THREAD_LIB});   
else
    textcommands = strjoin({'mex -silent -largeArrayDims -outdir', OUTDIR, ADD_INC, CUDA_INC...
    , mfile, strjoin(varargin), ['-L' path ' -lfftw3f'], ['-L' path ' -lfftw3'], FFTW_THREAD_LIB}); 
end;
disp(textcommands);
eval(textcommands);