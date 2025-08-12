%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-17 23:06:05
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 20:58:43
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function startup(root_dir)
if nargin < 1
    root_dir = '.';
end

warning('off','all')
addpath(genpath([root_dir, '/utils']))
addpath(genpath([root_dir,'/mexfun']))
% addpath(genpath([root_dir,'/sdpt3fun']))
addpath(genpath([root_dir,'/ALPDSN']))
addpath(genpath([root_dir,'/base']))
% addpath([root_dir,'linalg'])

% addpath_suitesparse;

warning('on','all')