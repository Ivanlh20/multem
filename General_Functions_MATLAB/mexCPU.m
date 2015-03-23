function [] = mexCPU(opt, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
  varargin{k} = strcat(path, '/', char(varargin{k}));
end;
CUDA_PATH = getenv('CUDA_PATH');
if(isempty(CUDA_PATH))
  CUDA_INC='-I/opt/cuda/include';
else
  CUDA_INC = strcat('-I"', CUDA_PATH, filesep, 'include"');
end;
if (opt==1)
    textcommands = strjoin({'mex -silent -g -largeArrayDims -outdir ../', strcat('-I', path), CUDA_INC, mfile, strjoin(varargin), ['-L' path ' -lfftw3']});   
else
    textcommands = strjoin({'mex -silent -largeArrayDims -outdir ../', strcat('-I', path), CUDA_INC, mfile, strjoin(varargin), ['-L' path ' -lfftw3']});
end;
disp(textcommands);
eval(textcommands);
