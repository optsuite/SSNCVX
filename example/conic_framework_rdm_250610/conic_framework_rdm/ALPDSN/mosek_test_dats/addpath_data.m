%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-13 16:15:57
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-21 09:54:44
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function data_dir = addpath_data(datasetname)
    if nargin < 1
        datasetname = '';
    else
        datasetname = [char(datasetname), '/'];
    end
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        % data_dir = ['/Users/niehantao/Desktop/software/data/sdp_data/', datasetname];
        data_dir = ['../../datasets', datasetname];
    elseif strcmp(computer, 'GLNXA64') % my linux server
        data_dir = ['/bicmr/home/optsuite/wt/datasets', datasetname];
    elseif strcmp(computer,'PCWIN64')
        computerName = getenv('COMPUTERNAME');
        if strcmp(computerName, "reserved for Zhanwang Deng")
            data_dir = ['E:\seafile\Seafile\ALDSDP\code\ssnsdp-beta-2nd\sdp_data\',datasetname];
        elseif strcmp(computerName, 'MAGICBOOK-TAL')
            data_dir = ['D:\Dev\OptSuite\datasets', datasetname];
        elseif strcmp(computerName, 'PC2024-TAL')
            data_dir = ['C:\Dev\optsuite\datasets', datasetname];
        else
            data_dir = ['../../datasets', datasetname];
        end
    end
end