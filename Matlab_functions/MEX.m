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

ADD_INC = strcat('-I', path);
OUTDIR = strcat('..', filesep, 'Mex_Executables');

if (strcmpi(option, 'debug'))
    textcommands = strjoin({'mex -silent -g -largeArrayDims -outdir', OUTDIR, ADD_INC, CUDA_INC...
    , mfile, strjoin(varargin), ['-L' path ' -llibfftw3-3'], ['-L' path ' -llibfftw3f-3']});   
else
    textcommands = strjoin({'mex -silent -largeArrayDims -outdir', OUTDIR, ADD_INC, CUDA_INC...
    , mfile, strjoin(varargin), ['-L' path ' -llibfftw3-3'], ['-L' path ' -llibfftw3f-3']}); 
end;
disp(textcommands);
eval(textcommands);