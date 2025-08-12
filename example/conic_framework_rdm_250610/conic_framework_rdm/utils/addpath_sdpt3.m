%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-20 17:37:16
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function addpath_sdpt3
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64')  % my mac laptop
        SDPT3Home = '/Users/niehantao/Desktop/software/SDPT3-4.0';
        eval(['addpath ',strcat(SDPT3Home,'/')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver/Mexfun')]);
        eval(['addpath ',strcat(SDPT3Home,'/HSDSolver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Examples')]);
    elseif strcmp(computer, 'GLNXA64') % my linux server
        SDPT3Home = '/bicmr/home/nieht/software/SDPT3-4.0';
        eval(['addpath ',strcat(SDPT3Home,'/')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver/Mexfun')]);
        eval(['addpath ',strcat(SDPT3Home,'/HSDSolver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Examples')]);
    end
end
