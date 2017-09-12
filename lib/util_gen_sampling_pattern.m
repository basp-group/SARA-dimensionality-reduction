function [uw, vw, Nm] = util_gen_sampling_pattern(pattern, param)

if ~isfield(param, 'fpartition') && ~isfield(param, 'fpartition_x') && ~isfield(param, 'fpartition_y')
    % symetric partitioning
    param.fpartition = [-pi pi];
end
if isfield(param, 'fpartition')
    % symetric partitioning
    param.fpartition_y = param.fpartition;
    param.fpartition_x = param.fpartition;
end


if strcmp(pattern, 'gaussian')
    sigma_m = param.g_sigma;
    Nm = round(param.p * param.N);
    
    u = sigma_m * randn(Nm, 1);
    v = sigma_m * randn(Nm, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    while length(sf) < Nm
        Nmextra = 2 * (Nm - length(sf));
        u = [u; sigma_m * randn(Nmextra, 1)];
        v = [v; sigma_m * randn(Nmextra, 1)];
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((u<pi) & (u>-pi));
        sfv = find((v<pi) & (v>-pi));
        sf = intersect(sfu, sfv);
    end
    
    vw = v(sf(1:Nm));
    uw = u(sf(1:Nm));
    
    
    Nm = length(uw);
end

if strcmp(pattern, 'gaussian-x2')
    sigma_m_1 = param.gx2_sigma_1;
    sigma_m_2 = param.gx2_sigma_2;
    Nm_1 = round(param.p * param.N * param.gx2_ratio);
    Nm_2 = round(param.p * param.N) - Nm_1;
    
    
    u = sigma_m_1 * randn(Nm_1, 1);
    v = sigma_m_1 * randn(Nm_1, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    while length(sf) < Nm_1
        Nmextra = 2 * (Nm_1 - length(sf));
        u = [u; sigma_m_1 * randn(Nmextra, 1)];
        v = [v; sigma_m_1 * randn(Nmextra, 1)];
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((u<pi) & (u>-pi));
        sfv = find((v<pi) & (v>-pi));
        sf = intersect(sfu, sfv);
    end
    
    vw = v(sf(1:Nm_1));
    uw = u(sf(1:Nm_1));
    
    
    
    u = sigma_m_2 * randn(Nm_2, 1);
    v = sigma_m_2 * randn(Nm_2, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    while length(sf) < Nm_2
        Nmextra = 2 * (Nm_2 - length(sf));
        u = [u; sigma_m_2 * randn(Nmextra, 1)];
        v = [v; sigma_m_2 * randn(Nmextra, 1)];
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((u<pi) & (u>-pi));
        sfv = find((v<pi) & (v>-pi));
        sf = intersect(sfu, sfv);
    end
    
    vw = [vw; v(sf(1:Nm_2))];
    uw = [uw; u(sf(1:Nm_2))];
    
    
    Nm = length(uw);
end

if strcmp(pattern, 'general-gaussian')
    Nm = round(param.p * param.N);
    
    uv = util_rand_2dggd(param.ggd_rho, param.ggd_beta, Nm);
    u = uv(:, 1);
    v = uv(:, 2);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    while length(sf) < Nm
        Nmextra = 2 * (Nm - length(sf));
        uv = util_rand_2dggd(param.ggd_rho, param.ggd_beta, Nmextra);
        u = [u; uv(:, 1)];
        v = [v; uv(:, 2)];
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((u<pi) & (u>-pi));
        sfv = find((v<pi) & (v>-pi));
        sf = intersect(sfu, sfv);
    end
    
    vw = v(sf(1:Nm));
    uw = u(sf(1:Nm));

    
    Nm = length(uw);
end

if strcmp(pattern, 'gaussian+large-holes')
    Nh = param.gh_hole_number;
    sigma_h = param.gh_sigma_holes;
    
    hu = [];
    hv = [];
    while length(hu) < Nh
        uv = -pi + 2*pi * rand(2, 1);
        % generate holes in the coverage
        if pdf('norm', 0, 0, sigma_h) * rand(1, 1) > pdf('norm', norm(uv), 0, sigma_h)
            hu = [hu; uv(1)];
            hv = [hv; uv(2)];
        end
    end

    
    % generate points outside the holes
    sigma_m = param.gh_sigma;
    Nm = round(param.p * param.N);
    
    u = sigma_m * randn(Nm, 1);
    v = sigma_m * randn(Nm, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    hs = param.gh_hole_size;
    for k = 1:Nh
        % discard points inside the holes
        sfu = find((u<hu(k)+hs) & (u>hu(k)-hs));
        sfv = find((v<hv(k)+hs) & (v>hv(k)-hs));
        sfh = intersect(sfu, sfv);
        sf = setdiff(sf, sfh);
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %d frequency points out of %d ... \n', length(u), Nm);
    
    while length(u) < Nm
        Nmextra = 2 * (Nm - length(u));
        us = sigma_m * randn(Nmextra, 1);
        vs = sigma_m * randn(Nmextra, 1);
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((us<pi) & (us>-pi));
        sfv = find((vs<pi) & (vs>-pi));
        sf = intersect(sfu, sfv);
        
        
        for k = 1:Nh
            % discard points inside the holes
            sfu = find((us<hu(k)+hs) & (us>hu(k)-hs));
            sfv = find((vs<hv(k)+hs) & (vs>hv(k)-hs));
            sfh = intersect(sfu, sfv);
            sf = setdiff(sf, sfh);
        end
        
        u = [u; us(sf)];
        v = [v; vs(sf)];
        fprintf('Computed %d frequency points out of %d ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %d frequency points out of %d ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    

    Nm = length(uw);
end

if strcmp(pattern, 'file')
    [uw, vw, ~, Nm] = util_uv_read(param.f_file_name);
end

if strcmp(pattern, 'file+undersample')
    [uw, vw, ~, Nm] = util_uv_read(param.fu_file_name);
    Nmn = round(param.p * param.N);
    
    if strcmp(param.fu_undersampling_type, 'uniform')    
        if Nmn > Nm
            error('Can''t undersample the UV coverage: Not enough points');
        end
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            Nm = Nm - length(r);
        end
    end
    
    if strcmp(param.fu_undersampling_type, 'gaussian')
        Nmo = length(uw);
        
        p_max = mvnpdf([0 0],[0 0], [param.fu_g_sigma 0; 0 param.fu_g_sigma]);
        prob = p_max * rand(Nmo, 1);

        sel = (prob > mvnpdf([uw vw],[0 0], [param.fu_g_sigma 0; 0 param.fu_g_sigma]));
        uw(sel) = [];
        vw(sel) = [];
        
        fprintf('Selected %d points out of %d using a Gaussian profile \n', length(uw), Nmo);
        
        Nm = length(uw);
        if Nmn > Nm
            error('Can''t undersample the UV coverage. \n Not enough points left after generating the Gaussian profile: only %d points left after resampling %d points', Nm, Nmn);
        end
        
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            Nm = Nm - length(r);
        end
        fprintf('Keeping only %d frequency points out of %d ... \n', Nmn, Nmo);
    end
    
    if strcmp(param.fu_undersampling_type, 'general-gaussian')
        Nmo = length(uw);
        mggdpdf = @(x, M)  exp(-(sum((x * diag(1./M)) .* x, 2)).^param.fu_ggd_beta/(2*param.fu_ggd_rho^param.fu_ggd_beta));
        
        p_max = mggdpdf([0 0], [param.fu_g_sigma param.fu_g_sigma]);
        prob = p_max * rand(Nmo, 1);

        sel = (prob > mggdpdf([uw vw], [param.fu_g_sigma param.fu_g_sigma]));
        uw(sel) = [];
        vw(sel) = [];
        
        fprintf('Selected %d points out of %d using a generalised Gaussian profile \n', length(uw), Nmo);
        
        Nm = length(uw);
        if Nmn > Nm
            error('Can''t undersample the UV coverage. \n Not enough points left after generating the Gaussian profile: only %d points left after resampling %d points', Nm, Nmn);
        end
        
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            Nm = Nm - length(r);
        end
        fprintf('Keeping only %d frequency points out of %d ... \n', Nmn, Nmo);
    end
end

if strcmp(pattern, 'gaussian+missing-pixels')
    fprintf('\n**** WARNING **** \nThis can be very slow.\n');
    prob_th = param.gmp_hole_prob;
    Noy = param.Noy;
    Nox = param.Nox;
    
    [hu, hv] = meshgrid(1:Noy, 1:Nox);
    hu = hu(:);
    hv = hv(:);
    
    prob = rand(param.Noy * param.Nox, 1);
    hu = hu(prob < prob_th);
    hv = hv(prob < prob_th);
    Nh = length(hu);
    
    
    hup = (2*pi*hu/Noy - pi) + pi/Noy;
    hum = (2*pi*hu/Noy - pi) - pi/Noy;
    
    hvp = (2*pi*hv/Nox - pi) + pi/Nox;
    hvm = (2*pi*hv/Nox - pi) - pi/Nox;
    
    % generate points outside the holes
    sigma_m = param.gmp_sigma;
    Nm = round(param.p * param.N);
    
    u = sigma_m * randn(Nm, 1);
    v = sigma_m * randn(Nm, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    for k = 1:Nh
        % discard points inside the holes        
        sfu = find(u<hup(k) & u>hum(k));
        sfv = find(v<hvp(k) & v>hvm(k));
        sfh = intersect(sfu, sfv);
        sf = setdiff(sf, sfh);
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %d frequency points out of %d ... \n', length(u), Nm);
    
    while length(u) < Nm
        Nmextra = 2 * (Nm - length(u));
        us = sigma_m * randn(Nmextra, 1);
        vs = sigma_m * randn(Nmextra, 1);
        
        
        % discard points outside (-pi,pi)x(-pi,pi)
        sfu = find((us<pi) & (us>-pi));
        sfv = find((vs<pi) & (vs>-pi));
        sf = intersect(sfu, sfv);
        
        
        for k = 1:Nh
            % discard points inside the holes        
            sfu = find(u<hup(k) & u>hum(k));
            sfv = find(v<hvp(k) & v>hvm(k));
            sfh = intersect(sfu, sfv);
            sf = setdiff(sf, sfh);        
        end
        
        u = [u; us(sf)];
        v = [v; vs(sf)];
        fprintf('Computed %d frequency points out of %d ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %d frequency points out of %d ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    

    Nm = length(uw);
end


end