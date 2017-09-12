function [epsilont, epsilonts, epsilon, epsilons] = util_gen_L2_bounds(y0, input_snr, sigma_noise, l2_ball_definition, stopping_criterion, equal_bounds, param)
% generates the input data


R = length(y0);
Nm = numel(cell2mat(y0));


if ~exist('equal_bounds', 'var')
    equal_bounds = 0;
end

if length(param.val_eps_v) == 1
    param.val_eps_v = ones(R, 1) * param.val_eps_v / sqrt(R);
end

if isempty(sigma_noise)
    % Gaussian i.i.d. noise
    normy0 = norm(cell2mat(y0));
    sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);
    fprintf('\nComputing sigma_noise from input SNR ... \n');
    fprintf('\nsigma_noise = %f \n', sigma_noise);
else
    fprintf('\nUsing provided sigma_noise ... \n');
    fprintf('\nsigma_noise = %f \n', sigma_noise);
end


if strcmp(l2_ball_definition, 'value')
    % estimate L2 ball parameter
    epsilon = norm(param.val_eps_v);
    epsilont = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater than the mean by
        % ~sqrt(mean(length(y{:})) (if equal length)
        epsilont{q} = param.val_eps_v(q);
    end
end

if strcmp(l2_ball_definition, 'sigma')
    s1 = param.sigma_ball;
    % estimate L2 ball parameter
    epsilon = sqrt(Nm + s1*sqrt(2*Nm)) * sigma_noise;
    epsilont = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater than the mean by
        % ~sqrt(mean(length(y{:})) (if equal length)
        epsilont{q} = sqrt(size(y0{q}, 1) + s1*sqrt(2*size(y0{q}, 1))) * sigma_noise;
    end
end

if strcmp(stopping_criterion, 'sigma')
    s2 = param.sigma_stop;
    % estimate L2 ball parameter
    epsilons = sqrt(Nm + s2*sqrt(2*Nm)) * sigma_noise;
    epsilonts = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater than the mean by
        % ~sqrt(mean(length(y{:})) (if equal length)
        epsilonts{q} = sqrt(size(y0{q}, 1) + s2*sqrt(2*size(y0{q}, 1))) * sigma_noise;
    end
end

if strcmp(l2_ball_definition, 'chi-percentile')
    p1 = param.chi_percentile_ball;
    % estimate L2 ball parameter
    epsilon = sqrt(chi2inv(p1, Nm)) * sigma_noise;
    epsilont = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater
        epsilont{q} = sqrt(chi2inv(p1, size(y0{q}, 1))) * sigma_noise;
    end
end

if strcmp(stopping_criterion, 'chi-percentile')
    p2 = param.chi_percentile_stop;
    % estimate L2 ball parameter
    epsilons = sqrt(chi2inv(p2, Nm)) * sigma_noise;
    epsilonts = cell(R,1);
    for q = 1:R
       % this produces a global bound which is greater
        epsilonts{q} = sqrt(chi2inv(p2, size(y0{q}, 1))) * sigma_noise;
    end
end

if strcmp(stopping_criterion, 'l2-ball-percentage' )
    % estimate L2 ball parameter
    sp = param.l2_ball_percentage_stop;
    epsilons = epsilon * sp;
    epsilonts = cell(R,1);
    for q = 1:R
       % this produces a global bound which is greater
        epsilonts{q} = epsilont{q} * sp;
    end
end


if equal_bounds
    n1 = norm(cell2mat(epsilont))/epsilon;
    n2 = norm(cell2mat(epsilonts))/epsilons;
    for q = 1:R
        epsilont{q} = epsilont{q} / n1;
    end

    for q = 1:R
        epsilonts{q} = epsilonts{q} / n2;
    end
end

end

