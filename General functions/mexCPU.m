function [] = mexCPU(opt, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
  varargin{k} = strcat(path, char(varargin{k}));
end;

if (opt==1)
    textcomands = strjoin({'mex -g', mfile, strjoin(varargin)});      
else
    textcomands = strjoin({'mex ', mfile, strjoin(varargin)});
end;

eval(textcomands);