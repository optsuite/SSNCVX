function [init] = init_parameters(params)
init = struct();
init.tau1 = 1;
init.tau2 = 1;
init.tau3 = 1;
init.tau4 = 1;
init.taux1 = 1;
init.taux2 = 1;
init.taux3 = 1;
init.taux4 = 1;

if isfield(params,'f') &&~isempty(params.f) && isempty(params.C)&& isempty(params.Q)

% elseif isscalar(params.pblk{1}) && strcmp(params.pblk{1}.type,'l1') && strcmp(params.f{1}.type,'l1')  %% Lasso
%     init.tau1 = 1;
%     init.tau2 = 1;
%     init.tau3 = 1;
%     init.tau4 = 1;
%     init.taux1 = 1;
%     init.taux2 = 1;
%     init.taux3 = 1;
%     init.taux4 = 1;
end





end