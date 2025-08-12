function [stepsize, smooth_mu] = linesearch(X, y, S, dx, dy, model, trans, param_ls)

    if nargin < 8
        param_ls = 1;
    end

    C = model.C;
    stepsize = trans.stepsize;
    best_stepsize = stepsize;
    energy_old = compute_dinf(y, S, trans);
    energy_best = inf;
    stepsize_min = trans.linesearch_minstep;
    rho = trans.linesearch_rho;
    dinf_threshold = trans.linesearch_threshold_dinf;
    ls_const = trans.ls_const;

    if dinf_threshold > 0 && energy_old > dinf_threshold
        smooth_mu = trans.smooth_mu;
        return;
    end

    W0 = trans.ATmap(y.var) - (C - X.var / trans.sigma);
    tmp = trans.ATmap(dy) + dx / trans.sigma;

    while true
        ynew = struct();
        ynew.var = y.var + dy * stepsize;
        
        W = W0 + tmp * stepsize;
        mu = 0;
        if trans.smooth
            mu = trans.smooth_mu / trans.sigma;
        elseif trans.smooth_linesearch_update_mu
            mu = trans.smooth_mu * trans.smooth_linesearch_mu_fact / trans.sigma;
        end

        if mu == 0
            [S0] = projectmit2(model.K, W);
        else
            [S0] = projectmit2smooth(model.K, W, mu, trans.socp_formula);
        end

        S.var = S0 - W;

        energy_temp = compute_dinf(ynew, S, trans);
        % fprintf('stepsize: %f, energy_temp: %f\n', stepsize, energy_temp);
        if energy_temp < energy_best
            energy_best = energy_temp;
            best_stepsize = stepsize;
        end
        if energy_best<= param_ls * ls_const * energy_old|| stepsize <= stepsize_min
            break;
        else
            stepsize= stepsize * rho;
        end
    end
    stepsize = best_stepsize;
    if ~trans.smooth
        smooth_mu = 0;
    elseif trans.smooth_linesearch_update_mu
        smooth_mu = trans.smooth_mu * trans.smooth_linesearch_mu_fact;
    else
        smooth_mu = trans.smooth_mu;
    end
end

function dinforg = compute_dinf(y, S, trans)
    y = y.var;
    y = trans.scale.Cscale * (trans.scale.DA * bwsolve(trans.Lchol, y));
    AtymCSnrm = norm((Atyfun(trans.K, trans.Atorg, y) + trans.scale.Cscale * S.var) - trans.Corg);
    dinforg = AtymCSnrm / (1 + norm(trans.Corg));
end

function pinforg = compute_pinf(X, trans)
    if ~isfield(X, 'Avarorg')
        [~, AXorg] = trans.Amap(X.var);
    else
        AXorg = X.Avarorg;
    end
    bmAXnrm = norm(AXorg ./ diag(trans.scale.DA) * trans.scale.bscale - trans.borg);
    pinforg = bmAXnrm / (1 + norm(trans.borg));
end
