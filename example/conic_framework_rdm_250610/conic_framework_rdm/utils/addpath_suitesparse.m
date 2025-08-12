%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-20 17:37:31
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function addpath_suitesparse
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        run("/Users/niehantao/Desktop/software/SuiteSparse/SuiteSparse_paths")
    elseif strcmp(computer, 'GLNXA64') % my linux server
        run("/bicmr/home/nieht/software/SuiteSparse/SuiteSparse_paths");
    end
end

