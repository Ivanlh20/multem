function [] = mexGPU(opt, mfile, path, varargin)
nVarargs = length(varargin);
for k = 1:nVarargs
  varargin{k} = strcat(path, char(varargin{k}));
end;

ext = '"-IC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v6.5\include" "-LC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v6.5\lib\x64" -lcuda -lcudart -lcufft';
if (opt==1)
    textcomands = strjoin({'nvmex -g -f nvmexopts64.bat', mfile, strjoin(varargin), ext});      
else
    textcomands = strjoin({'nvmex -f nvmexopts64.bat', mfile, strjoin(varargin), ext});
end;
eval(textcomands);

% function [] = mexGPU(opt, mfile, path, varargin)
% nVarargs = length(varargin);
% for k = 1:nVarargs
%   varargin{k} = strcat(path, char(varargin{k}));
% end;
% 
% ext = '-IC\cuda\include -LC:\cuda\lib\x64 -lcuda -lcudart -lcufft';
% if (opt==1)
%     textcomands = strjoin({'nvmex -g -f nvmexopts64.bat', mfile, strjoin(varargin), ext});      
% else
%     textcomands = strjoin({'nvmex -f nvmexopts64.bat', mfile, strjoin(varargin), ext});
% end;
% eval(textcomands);