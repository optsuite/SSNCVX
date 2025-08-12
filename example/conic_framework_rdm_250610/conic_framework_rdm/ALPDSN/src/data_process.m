function [bnew,At,cnz] = data_process(blk,At,b)
% DATA_PROCESS: a data process code.
%
% Copyright (c) 2017 by
% Yongfeng Li, Zaiwen Wen
%%



b(abs(b)<1e-13) = 0;

nblk = size(blk,1);
cnz = zeros(size(b));
for i = 1:nblk
    pblk = blk(i,:);
    if length(pblk{2}) ==  1 && strcmp(pblk{1},'s')
        spa = mexsmat_sdpnal(pblk,At{i}*b);
        pm = symamd(spa);
        Z = spa(pm,pm);
        if isdiag(Z)
            break;
        end

        [t,~] = etree(Z);
        idx0 = find(t == 0);
        len0 = length(idx0);
        nzpm = zeros(size(spa));
        

        offset = 1;
        for it = 1:len0
            idx = offset:idx0(it);
            nit = length(idx);
            matpm = pm(idx);
            matpm = sort(matpm,'ascend');
            nzpm(matpm,matpm) = randn(nit,nit);
            offset = idx0(it)+1;
        end

        bb = At{i}'*mexsvec_sdpnal(pblk,nzpm);
        cnz = (bb ~= 0) | cnz;
    end
end

cnz0 = b ~= 0;
cnz = cnz0 | cnz;
bnew = b(cnz);
for p = 1:length(At)
    At{p} = At{p}(:,cnz);
end
