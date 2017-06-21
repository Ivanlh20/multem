function [] = mexBlasLapackCPU(opt, mfile, path, varargin)
    function [so] = qt(si)
        so = strcat('''', si, '''');
    end

    nVarargs = length(varargin);
    for k = 1:nVarargs
      varargin{k} = strcat('''', strcat(path,'/'), char(varargin{k}), '''');
    end

    lapacklib = fullfile(matlabroot, ...
      'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
    lapacklib = qt(lapacklib);

    blaslib = fullfile(matlabroot, ...
      'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
    blaslib = qt(blaslib);
    
    CUDA_INC = strcat('-I"', getenv('CUDA_PATH'), filesep, 'include"');
    mfile = qt(mfile);
    if (opt==1)
        textcomands = strjoin({qt('-g'), qt('-largeArrayDims'), qt(strcat('-I', path)), qt(CUDA_INC), mfile, strjoin(varargin, ','), blaslib, lapacklib},',');     
    else
        textcomands = strjoin({qt('-largeArrayDims'), qt(strcat('-I', path)), qt(CUDA_INC), mfile, strjoin(varargin, ','), blaslib, lapacklib},',');       
    end
    textcomands = strcat('mex(', textcomands, ')');
    eval(textcomands);
end


% lapacklib = fullfile(matlabroot, ...
%   'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
% blaslib = fullfile(matlabroot, ...
%   'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');