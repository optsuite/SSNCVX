%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-12 12:15:41
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function addpath_mosek
    if strcmp(computer, "MACI64")
        addpath(genpath('/Users/niehantao/Desktop/software/mosek_maci/10.0/toolbox/r2017a'))
    elseif strcmp(computer, "MACA64")
            addpath(genpath('/Users/niehantao/Desktop/software/mosek_maca/10.1/toolbox/r2022b'))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath("/bicmr/home/nieht/software/mosek/10.0/toolbox/r2017a");
    end
end