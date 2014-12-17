function [] = mexCPU(opt, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
  varargin{k} = strcat(path, '/', char(varargin{k}));
end;
CUDA_INC = strcat('-I"', getenv('CUDA_PATH'), filesep, 'include"');
if (opt==1)
    textcomands = strjoin({'mex -g -largeArrayDims', strcat('-I', path), CUDA_INC, '-I/opt/cuda/include', mfile, strjoin(varargin)});   
else
    textcomands = strjoin({'mex -largeArrayDims', strcat('-I', path), CUDA_INC, '-I/opt/cuda/include', mfile, strjoin(varargin)});
end;

eval(textcomands);