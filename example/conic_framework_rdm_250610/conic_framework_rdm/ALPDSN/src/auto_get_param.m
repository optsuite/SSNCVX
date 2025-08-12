function [opts] = auto_get_param(opts, m, sp_info)
    if isempty(sp_info.row_sp) && isempty(sp_info.col_sp)
        disp('组合一:无稀疏行列');
        opts.linesearch = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1e-06;
        opts.smooth_threshold = 1e-13;
        opts.resratio = 0.1;
        opts.smooth_ratio_update_mu = 1;
    
        opts.muopts.mu_update_itr = 1000;
        opts.gamma1 = 0.01;
    
        opts.max_non_monotone = 20;
        opts.param_am = 0.8;
        opts.adaptive_mu = 1;
        opts.smooth_linesearch_mu_fact = 0.5;
        
    elseif isempty(sp_info.row_den) && isempty(sp_info.col_den)
        disp('组合二:无稠密行列');
        opts.smooth_mu = 1e-05;
        opts.gamma1 = 0.01;
        opts.resratio = 1;
        opts.smooth_ratio_update_mu = 0.5;
        
    elseif isempty(sp_info.col_den)
        disp('组合三:有稠密行无稠密列');
        opts.linesearch = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1e-06;
        opts.smooth_threshold = 1e-15;
        opts.resratio = 1;
        opts.smooth_ratio_update_mu = 1;
        opts.linesearch_min_step = 0.1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
    
        opts.muopts.mu_update_itr = 1000;
        opts.gamma1 = 0.3;
        opts.gamma3 = 1;
    
        opts.adaptive_mu = 1;
        opts.smooth_linesearch_mu_fact = 0.7;
        opts.max_non_monotone = 8;
        opts.param_am = 100;
        
    elseif ~isempty(sp_info.row_den) && ~isempty(sp_info.col_den)
        disp('组合四:有稠密行有稠密列');

        opts.scale_A = 0;
        opts.scale_bc_flag = 1;
        opts.adaptive_mu = 1;
        opts.gamma1 = 0.3;
        opts.gamma2 = 0.7;
        opts.gamma3 = 1;
        opts.linesearch = 1;
        opts.linesearch_min_step = 1.00e-01;
        opts.max_non_monotone = 4;
        opts.muopts.mu_fact = 0.6;
        opts.muopts.init = 100;
        opts.muopts.mu_update_itr = 10000;
        opts.muopts.mu_fact = 0.5;
        opts.param_am = 100;
        opts.resratio = 1;
        opts.sigxl = 0.3;
        opts.sigxu = 0.4;
        opts.sigyl = 0.6;
        opts.sigyu = 1;
        opts.smooth = 1;
        opts.smooth_linesearch_update_mu = 1;
        opts.smooth_mu = 1.00e-05;
        opts.smooth_threshold = 1.00e-15;
        opts.socp_formula = 3;
        

        % opts.linesearch = 1;
        % opts.smooth_linesearch_update_mu = 1;
        % opts.smooth_mu = 1e-05;
        % opts.smooth_threshold = 1e-15;
        % opts.resratio = 1;
        % opts.smooth_ratio_update_mu = 1;
        % opts.linesearch_min_step = 0.1;
        % opts.sigxl = 0.3;
        % opts.sigxu = 0.4;
        % opts.sigyl = 0.6;
        % opts.sigyu = 1;
        % opts.sigma = 10;
        % opts.scale_A = 0;
        % opts.scale_bc_flag = 0;
        % 
        % opts.muopts.mu_update_itr = 1000;
        % opts.gamma1 = 0.3;
        % opts.gamma3 = 1;
        % 
        % opts.adaptive_mu = 1;
        % opts.smooth_linesearch_mu_fact = 0.5;
        % opts.max_non_monotone = 50;
        % opts.param_am = 100;
        % opts.refinement_tol = 8;
        % opts.refine_max_iter = 3;
        
    elseif isempty(sp_info.row_den)
        if m < 1e+05
            disp('组合五:有稠密列无稠密行,m<1e+05');
            opts.linesearch = 1;
            opts.smooth_linesearch_update_mu = 1;
            opts.smooth_mu = 1e-01;
            opts.smooth_threshold = 1e-10;
            opts.resratio = 1;
            opts.smooth_ratio_update_mu = 1;
            opts.linesearch_min_step = 0.1;
            opts.sigxl = 0.3;
            opts.sigxu = 0.4;
            opts.sigyl = 0.6;
            opts.sigyu = 1;
            opts.sigma = 0.1;
    
            opts.muopts.mu_update_itr = 1000;
            opts.gamma1 = 0.3;
            opts.gamma3 = 1.5;
    
            opts.adaptive_mu = 1;
            opts.smooth_linesearch_mu_fact = 0.7;
            opts.max_non_monotone = 50;
            opts.param_am = 100;
            opts.refinement_tol = 13;
            opts.refine_max_iter = 1;
    
        elseif m < 2e+05
            disp('组合六:有稠密列无稠密行,1e+05<m<2e+05');
            opts.linesearch = 1;
            opts.smooth_linesearch_update_mu = 1;
            opts.smooth_mu = 4e-06;
            opts.smooth_threshold = 1e-16;
            opts.resratio = 1;
            opts.smooth_ratio_update_mu = 1;
            opts.linesearch_min_step = 1e-06;
            opts.sigma = 1e-03;
    
            opts.muopts.mu_update_itr = 1000;
            opts.gamma1 = 0.3;
    
            opts.adaptive_mu = 1;
            opts.smooth_linesearch_mu_fact = 0.3;
            opts.max_non_monotone = 5;
            opts.param_am = 10;
            opts.reorder_type = 3;
    
        elseif m < 1e+06
            disp('组合七:有稠密列无稠密行,2e+05<m<1e+06');
            opts.linesearch = 1;
            opts.smooth_linesearch_update_mu = 1;
            opts.smooth_mu = 1e-03;
            opts.smooth_threshold = 1e-14;
            opts.resratio = 0.5;
            opts.smooth_ratio_update_mu = 1;
            opts.linesearch_min_step = 1e-02;
            opts.sigyl = 0.5;
            opts.sigyu = 0.5;
            opts.sigma = 1;

            opts.scale_A = 0;
            opts.scale_bc_flag = 0;
    
            opts.muopts.mu_update_itr = 1000;
            opts.gamma1 = 0.3;
    
            opts.adaptive_mu = 1;
            opts.smooth_linesearch_mu_fact = 0.6;
            opts.max_non_monotone = 15;
            opts.param_am = 1000;
            opts.refinement_tol = 8;
            opts.refine_max_iter = 1;
    
        else
            disp('组合八:');
            opts.linesearch = 1;
            opts.smooth_linesearch_update_mu = 1;
            opts.smooth_mu = 1e-06;
            opts.smooth_threshold = 4e-16;
            opts.resratio = 0.5;
            opts.smooth_ratio_update_mu = 1;
            opts.linesearch_min_step = 0.1;
            opts.sigyl = 0.5;
            opts.sigyu = 0.5;
            opts.sigma = 3;
    
            opts.muopts.mu_update_itr = 1000;
            opts.gamma1 = 0.3;
    
            opts.adaptive_mu = 1;
            opts.smooth_linesearch_mu_fact = 0.5;
            opts.max_non_monotone = 5;
            opts.param_am = 1000;
            opts.refinement_tol = 13;
            opts.refine_max_iter = 3;
            opts.activate_nm = 1e-04;
            opts.non_monotone = 1;
            opts.param_nm = 5;
        end
    end
end

% function [trans, muopts, NEWTopts] = auto_get_param(trans, muopts, NEWTopts, m, sp_info)
%     if isempty(sp_info.row_sp) && isempty(sp_info.col_sp)
%         disp('组合一:无稀疏行列');
%         trans.linesearch = 1;
%         trans.smooth_linesearch_update_mu = 1;
%         trans.smooth_mu = 1e-06;
%         trans.smooth_threshold = 1e-13;
%         trans.resratio = 0.1;
%         trans.smooth_ratio_update_mu = 1;
% 
%         muopts.NEWT.mu_update_itr = 1000;
%         trans.muopts.NEWT.mu_update_itr = 1000;
%         NEWTopts.gamma1 = 0.01;
%         trans.NEWTopts.gamma1 = 0.01;
% 
%         trans.max_non_monotone = 20;
%         trans.param_am = 0.8;
%         trans.adaptive_mu = 1;
%         trans.smooth_linesearch_mu_fact = 0.5;
% 
%     elseif isempty(sp_info.row_den) && isempty(sp_info.col_den)
%         disp('组合二:无稠密行列');
%         trans.smooth_mu = 1e-05;
%         NEWTopts.gamma1 = 0.01;
%         trans.NEWTopts.gamma1 = 0.01;
%         trans.resratio = 1;
%         trans.smooth_ratio_update_mu = 0.5;
% 
%     elseif isempty(sp_info.col_den)
%         disp('组合三:有稠密行无稠密列');
%         trans.linesearch = 1;
%         trans.smooth_linesearch_update_mu = 1;
%         trans.smooth_mu = 1e-06;
%         trans.smooth_threshold = 1e-15;
%         trans.resratio = 1;
%         trans.smooth_ratio_update_mu = 1;
%         trans.linesearch_min_step = 0.1;
%         trans.sigxl = 0.3;
%         trans.sigxu = 0.4;
%         trans.sigyl = 0.6;
%         trans.sigyu = 1;
% 
%         muopts.NEWT.mu_update_itr = 1000;
%         trans.muopts.NEWT.mu_update_itr = 1000;
%         NEWTopts.gamma1 = 0.3;
%         trans.NEWTopts.gamma1 = 0.3;
%         NEWTopts.gamma3 = 1;
%         trans.NEWTopts.gamma3 = 1;
% 
%         trans.adaptive_mu = 1;
%         trans.smooth_linesearch_mu_fact = 0.7;
%         trans.max_non_monotone = 8;
%         trans.param_am = 100;
% 
%     elseif ~isempty(sp_info.row_den) && ~isempty(sp_info.col_den)
%         disp('组合四:有稠密行有稠密列');
%         trans.linesearch = 1;
%         trans.smooth_linesearch_update_mu = 1;
%         trans.smooth_mu = 1e-05;
%         trans.smooth_threshold = 1e-15;
%         trans.resratio = 1;
%         trans.smooth_ratio_update_mu = 1;
%         trans.linesearch_min_step = 0.1;
%         trans.sigxl = 0.3;
%         trans.sigxu = 0.4;
%         trans.sigyl = 0.6;
%         trans.sigyu = 1;
%         trans.sigma = 10;
% 
%         muopts.NEWT.mu_update_itr = 1000;
%         trans.muopts.NEWT.mu_update_itr = 1000;
%         NEWTopts.gamma1 = 0.3;
%         trans.NEWTopts.gamma1 = 0.3;
%         NEWTopts.gamma3 = 1;
%         trans.NEWTopts.gamma3 = 1;
% 
%         trans.adaptive_mu = 1;
%         trans.smooth_linesearch_mu_fact = 0.5;
%         trans.max_non_monotone = 50;
%         trans.param_am = 100;
%         trans.refinement_tol = 8;
%         trans.refine_max_iter = 3;
% 
%     elseif isempty(sp_info.row_den)
%         if m < 1e+05
%             disp('组合五:有稠密列无稠密行,m<1e+05');
%             trans.linesearch = 1;
%             trans.smooth_linesearch_update_mu = 1;
%             trans.smooth_mu = 1e-01;
%             trans.smooth_threshold = 1e-10;
%             trans.resratio = 1;
%             trans.smooth_ratio_update_mu = 1;
%             trans.linesearch_min_step = 0.1;
%             trans.sigxl = 0.3;
%             trans.sigxu = 0.4;
%             trans.sigyl = 0.6;
%             trans.sigyu = 1;
%             trans.sigma = 0.1;
% 
%             muopts.NEWT.mu_update_itr = 1000;
%             trans.muopts.NEWT.mu_update_itr = 1000;
%             NEWTopts.gamma1 = 0.3;
%             trans.NEWTopts.gamma1 = 0.3;
%             NEWTopts.gamma3 = 1.5;
%             trans.NEWTopts.gamma3 = 1.5;
% 
%             trans.adaptive_mu = 1;
%             trans.smooth_linesearch_mu_fact = 0.7;
%             trans.max_non_monotone = 50;
%             trans.param_am = 100;
%             trans.refinement_tol = 13;
%             trans.refine_max_iter = 1;
% 
%         elseif m < 2e+05
%             disp('组合六:有稠密列无稠密行,1e+05<m<2e+05');
%             trans.linesearch = 1;
%             trans.smooth_linesearch_update_mu = 1;
%             trans.smooth_mu = 4e-06;
%             trans.smooth_threshold = 1e-16;
%             trans.resratio = 1;
%             trans.smooth_ratio_update_mu = 1;
%             trans.linesearch_min_step = 1e-06;
%             trans.sigma = 1e-03;
% 
%             muopts.NEWT.mu_update_itr = 1000;
%             trans.muopts.NEWT.mu_update_itr = 1000;
%             NEWTopts.gamma1 = 0.3;
%             trans.NEWTopts.gamma1 = 0.3;
% 
%             trans.adaptive_mu = 1;
%             trans.smooth_linesearch_mu_fact = 0.3;
%             trans.max_non_monotone = 5;
%             trans.param_am = 10;
%             trans.reorder_type = 3;
% 
%         elseif m < 1e+06
%             disp('组合七:有稠密列无稠密行,2e+05<m<1e+06');
%             trans.linesearch = 1;
%             trans.smooth_linesearch_update_mu = 1;
%             trans.smooth_mu = 1e-03;
%             trans.smooth_threshold = 1e-14;
%             trans.resratio = 0.5;
%             trans.smooth_ratio_update_mu = 1;
%             trans.linesearch_min_step = 1e-02;
%             trans.sigyl = 0.5;
%             trans.sigyu = 0.5;
%             trans.sigma = 1;
% 
%             trans.scale_A = 0;
%             trans.scale_bc_flag = 0;
% 
%             muopts.NEWT.mu_update_itr = 1000;
%             trans.muopts.NEWT.mu_update_itr = 1000;
%             NEWTopts.gamma1 = 0.3;
%             trans.NEWTopts.gamma1 = 0.3;
% 
%             trans.adaptive_mu = 1;
%             trans.smooth_linesearch_mu_fact = 0.6;
%             trans.max_non_monotone = 15;
%             trans.param_am = 1000;
%             trans.refinement_tol = 8;
%             trans.refine_max_iter = 1;
% 
%         else
%             disp('组合八:');
%             trans.linesearch = 1;
%             trans.smooth_linesearch_update_mu = 1;
%             trans.smooth_mu = 1e-06;
%             trans.smooth_threshold = 4e-16;
%             trans.resratio = 0.5;
%             trans.smooth_ratio_update_mu = 1;
%             trans.linesearch_min_step = 0.1;
%             trans.sigyl = 0.5;
%             trans.sigyu = 0.5;
%             trans.sigma = 3;
% 
%             muopts.NEWT.mu_update_itr = 1000;
%             trans.muopts.NEWT.mu_update_itr = 1000;
%             NEWTopts.gamma1 = 0.3;
%             trans.NEWTopts.gamma1 = 0.3;
% 
%             trans.adaptive_mu = 1;
%             trans.smooth_linesearch_mu_fact = 0.5;
%             trans.max_non_monotone = 5;
%             trans.param_am = 1000;
%             trans.refinement_tol = 13;
%             trans.refine_max_iter = 3;
%             trans.activate_nm = 1e-04;
%             trans.non_monotone = 1;
%             trans.param_nm = 5;
%         end
%     end
% end