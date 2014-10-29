function [] = mexBlasLapackCPU(opt, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
  varargin{k} = strcat('''', path, char(varargin{k}), '''');
end;

lapacklib = fullfile(matlabroot, ...
  'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
lapacklib = strcat('''', lapacklib, '''');

blaslib = fullfile(matlabroot, ...
  'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
blaslib = strcat('''', blaslib, '''');

mfile = strcat('''', mfile, '''');
mopt = '-g'; mopt = strcat('''', mopt, '''');
if (opt==1)
    textcomands = strjoin({mopt, mfile, strjoin(varargin, ','), blaslib, lapacklib}, ',');  
    textcomands = strcat('mex(', textcomands, ')');      
else
    textcomands = strjoin({mfile, strjoin(varargin, ','), blaslib, lapacklib}, ',');  
    textcomands = strcat('mex(', textcomands, ')'); 
end;

eval(textcomands);

% lapacklib = fullfile(matlabroot, ...
%   'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
% blaslib = fullfile(matlabroot, ...
%   'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');