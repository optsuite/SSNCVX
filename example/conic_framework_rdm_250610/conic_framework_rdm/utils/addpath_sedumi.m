%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-20 17:37:21
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function addpath_sedumi
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        addpath("/Users/niehantao/Desktop/software/sedumi");
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath("/bicmr/home/nieht/software/sedumi");
    end
end
