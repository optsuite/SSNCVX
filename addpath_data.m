
function data_dir = addpath_data(datasetname)
    if nargin < 1
        datasetname = '';
    else
        datasetname = [char(datasetname), '/'];
    end
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        data_dir = ['./', datasetname];
    elseif strcmp(computer, 'GLNXA64') % my linux server
        data_dir = ['./', datasetname];
    elseif strcmp(computer,'PCWIN64')
        data_dir = ['./',datasetname];
    end
    addpath(genpath(data_dir));  
end
