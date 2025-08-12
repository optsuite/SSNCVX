function opts = default_opts(opts)
%% ------------------------------------------------------------------
if ~isfield(opts, 'log_path'); opts.log_path = ""; end
if ~isfield(opts, 'record'); opts.record = 1; end
if ~isfield(opts, 'scale_data'); opts.scale_data = 1; end
if ~isfield(opts, 'pmu'); opts.pmu = 1; end
if ~isfield(opts, 'sstruct'); opts.sstruct = 1; end
if ~isfield(opts, 'maxits'); opts.maxits = 1000; end
if ~isfield(opts, 'tol'); opts.tol = 1e-6; end
if ~isfield(opts, 'tolscale'); opts.tolscale = 1; end

%% -------------------------------------------------------------------------
% options for fixed-point algorithms
if ~isfield(opts, 'ADMMopts'); opts.ADMMopts = struct; end
ADMMopts = opts.ADMMopts;
if ~isfield(ADMMopts, 'print_itr'); ADMMopts.print_itr = 20; end
if ~isfield(ADMMopts, 'rho'); ADMMopts.rho = 1.618; end
if ~isfield(ADMMopts, 'imaxit'); ADMMopts.imaxit = 20000; end
opts.ADMMopts = ADMMopts;

%% ------------------------------------------------------------------------
% options for the semismooth Newton algorithm
if ~isfield(opts, 'NEWTopts'); opts.NEWTopts = struct; end
NEWTopts = opts.NEWTopts;
if ~isfield(NEWTopts,'sigPowx');   NEWTopts.sigPowx = 1;       end
if ~isfield(NEWTopts,'sigPowy');   NEWTopts.sigPowy = 1;       end
if ~isfield(NEWTopts,'sigPowz');   NEWTopts.sigPowz = 1;       end
if ~isfield(NEWTopts,'sigPowq');   NEWTopts.sigPowq = 1;       end
if ~isfield(NEWTopts, 'resFac'); NEWTopts.resFac = 0.98; end
if ~isfield(NEWTopts, 'eta1'); NEWTopts.eta1 = 1e-4; end
if ~isfield(NEWTopts, 'eta2'); NEWTopts.eta2 = 0.9; end
if ~isfield(NEWTopts, 'gamma1'); NEWTopts.gamma1 = opts.gamma1; end %可调
if ~isfield(NEWTopts, 'gamma2'); NEWTopts.gamma2 = opts.gamma2; end %可调
if ~isfield(NEWTopts, 'gamma3'); NEWTopts.gamma3 = opts.gamma3; end %可调
if ~isfield(NEWTopts, 'lambda'); NEWTopts.lambda = 1; end %adjust mu
if ~isfield(NEWTopts, 'print_itr'); NEWTopts.print_itr = 1; end
opts.NEWTopts = NEWTopts;

%% ------------------------------------------------------------------------
% parameter for newton system
if ~isfield(opts, 'method'); opts.method = 'iterative'; end


%% ------------------------------------------------------------------------
% parameter for CG
if ~isfield(opts, 'cgopts'); opts.cgopts = struct; end
cgopts = opts.cgopts;
if ~isfield(cgopts, 'CG_maxit'); cgopts.CG_maxit = 500; end %可调
if ~isfield(cgopts, 'CG_tol'); cgopts.CG_tol = 1e-2; end
if ~isfield(cgopts, 'cgtolmin'); cgopts.cgtolmin = opts.cgtol; end
if ~isfield(cgopts, 'CG_adapt'); cgopts.CG_adapt = 1; end
opts.cgopts = cgopts;

%% ------------------------------------------------------------------
% parameters for adjusting pmu
if ~isfield(opts, 'muopts'); opts.muopts = struct; end
muopts = opts.muopts;
if ~isfield(muopts, 'adp_mu'); muopts.adp_mu = 1; end %1 or 0
if ~isfield(muopts, 'NEWT'); muopts.NEWT = struct; end
if ~isfield(muopts.NEWT, 'mu_min'); muopts.NEWT.mu_min = 1e-6; end
if ~isfield(muopts.NEWT, 'mu_max'); muopts.NEWT.mu_max = 1e3; end
if ~isfield(muopts.NEWT, 'mu_update_itr'); muopts.NEWT.mu_update_itr = muopts.mu_update_itr; end %10 可调
if ~isfield(muopts.NEWT, 'mu_delta'); muopts.NEWT.mu_delta = 5e-1; end % 可调
if ~isfield(muopts.NEWT, 'mu_fact'); muopts.NEWT.mu_fact = muopts.mu_fact; end %5/3 可调


opts.muopts = muopts;
end