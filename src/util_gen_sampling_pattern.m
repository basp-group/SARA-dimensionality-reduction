function [u, v, uw, vw, Nm] = util_gen_sampling_pattern(pattern, param)

if ~isfield(param, 'fpartition') && ~isfield(param, 'fpartition_x') && ~isfield(param, 'fpartition_y')
    % symetric partitioning
    param.fpartition = [-pi pi];
end
if isfield(param, 'fpartition')
    % symetric partitioning
    param.fpartition_y = param.fpartition;
    param.fpartition_x = param.fpartition;
end
if ~isfield(param, 'sigma_holes')
    % variance of the inverse gaussian for selecting the missing pixels
    param.sigma_holes = pi/2;
end
if ~isfield(param, 'hole_number')
    % filling level of the image
    param.hole_number = 10;
end
if ~isfield(param, 'hole_size')
    % filling level of the image
    param.hole_size = pi/8;
end
if ~isfield(param, 'sigma')
    % variance of the gaussian for generating measurements
    param.sigma = pi/3;
    param.fpartition_x = param.fpartition;
end

if strcmp(pattern, 'gaussian')
    sigma_m = param.sigma;
    Nm = ceil(param.p * param.N);
    
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
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end

if strcmp(pattern, 'gaussian+large-holes')
    Nh = param.hole_number;
    sigma_h = param.sigma_holes;
    
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
    sigma_m = param.sigma;
    Nm = param.p * param.N;
    
    u = sigma_m * randn(Nm, 1);
    v = sigma_m * randn(Nm, 1);

    % discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    hs = param.hole_size;
    for k = 1:Nh
        % discard points inside the holes
        sfu = find((u<hu(k)+hs) & (u>hu(k)-hs));
        sfv = find((v<hv(k)+hs) & (v>hv(k)-hs));
        sfh = intersect(sfu, sfv);
        sf = setdiff(sf, sfh);
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    
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
        fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %i frequency points out of %i ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end

if strcmp(pattern, 'uniform+large-holes')
    Nh = param.hole_number;
    sigma_h = param.sigma_holes;
    
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
    sigma_m = param.sigma;
    Nm = param.p * param.N;
    
    u = -pi + 2*pi * rand(Nm, 1);
    v = -pi + 2*pi * rand(Nm, 1);

% discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    hs = param.hole_size;
    for k = 1:Nh
        % discard points inside the holes
        sfu = find((u<hu(k)+hs) & (u>hu(k)-hs));
        sfv = find((v<hv(k)+hs) & (v>hv(k)-hs));
        sfh = intersect(sfu, sfv);
        sf = setdiff(sf, sfh);
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    
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
        fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %i frequency points out of %i ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end

if strcmp(pattern, 'geometric+large-holes')
    Nh = param.hole_number;
    sigma_h = param.sigma_holes;
    
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
    sigma_m = param.sigma;
    Nm = param.p * param.N;
    
    u = random('Geometric', 0.2, Nm, 1);
    v = random('Geometric', 0.2, Nm, 1);
    u = 2*pi* u/max(u); v = 2*pi* v/max(v);
    u = -pi + u; v = -pi + v;

% discard points outside (-pi,pi)x(-pi,pi)
    sfu = find((u<pi) & (u>-pi));
    sfv = find((v<pi) & (v>-pi));
    sf = intersect(sfu, sfv);
    
    hs = param.hole_size;
    for k = 1:Nh
        % discard points inside the holes
        sfu = find((u<hu(k)+hs) & (u>hu(k)-hs));
        sfv = find((v<hv(k)+hs) & (v>hv(k)-hs));
        sfh = intersect(sfu, sfv);
        sf = setdiff(sf, sfh);
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    
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
        fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %i frequency points out of %i ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end

if strcmp(pattern, 'file')
    [uw, vw, ~, Nm] = util_uv_read(param.file_name);
    
    u = uw;
    v = vw;
end

if strcmp(pattern, 'file+undersample')
    [uw, vw, ~, Nm] = util_uv_read(param.file_name);
    
    
    Nmn = round(param.p * param.N);
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
    u = uw;
    v = vw;
end



if strcmp(pattern, 'gaussian+missing-pixels')

    prob_th = param.hole_prob;
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
    sigma_m = param.sigma;
    Nm = param.p * param.N;
    
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
        if mod(k, 1000) == 0
            k
        end
    end
    
    u = u(sf);
    v = v(sf);    
    
    fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    
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

            if mod(k, 1000) == 0
                k
            end            
        end
        
        u = [u; us(sf)];
        v = [v; vs(sf)];
        fprintf('Computed %i frequency points out of %i ... \n', length(u), Nm);
    end
    
    fprintf('Keeping only %i frequency points out of %i ... \n', Nm, length(u));
    vw = v(1:Nm);
    uw = u(1:Nm);
    
    clear u;
    clear v;
    
    u = uw;
    v = vw;
    
    Nm = length(uw);
end



if isfield(param, 'fpartition_y') && isfield(param, 'fpartition_x')
    %% frequency clustering

    indexu = quantiz(uw, param.fpartition_y);
    indexv = quantiz(vw, param.fpartition_x);
    indexuv = length(param.fpartition_y)*indexu + indexv;
    indexuv = indexuv+1;

    hs = length(param.fpartition_x) * length(param.fpartition_y);

    u = cell(hs, 1);
    v = cell(hs, 1);
    for q = 1:hs
        u{q} = uw(indexuv==q);
        v{q} = vw(indexuv==q);
    end
end

if isfield(param, 'partition')
    %% data clustering

    hs = length(param.partition);

    u = cell(hs, 1);
    v = cell(hs, 1);
    Rp = 0;
    for q = 1:hs
        u{q} = uw(Rp+1:Rp+param.partition(q));
        v{q} = vw(Rp+1:Rp+param.partition(q));
        Rp = Rp + param.partition(q);
    end
end



end
