
function data_dir = addpath_data(datasetname)
    if nargin < 1
        datasetname = '';
    else
        datasetname = [char(datasetname), '/'];
    end
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        data_dir = ['/Users/niehantao/Desktop/software/data/sdp_data/', datasetname];
    elseif strcmp(computer, 'GLNXA64') % my linux server
        data_dir = ['/public/shared/sdp_data/', datasetname];
    elseif strcmp(computer,'PCWIN64')
        data_dir = ['../data/',datasetname];
end
