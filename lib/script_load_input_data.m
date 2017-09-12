% script that loads or generates the input data for the test

% set random nr gen
rng shuffle;

fprintf('Loading data from disk ... \n\n');
if exist(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)), 'file') || ... 
      exist(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)), 'file')

    vars = load(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)));
    y = vars.y;
    nWw = vars.nWw;
    yf = vars.yf;
    y0 = vars.y0;
    y0f = vars.y0f;

    vars = load(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)));
    
    num_tests = vars.num_tests;
    use_real_visibilities = vars.use_real_visibilities;
    input_snr = vars.input_snr;
    image_file_name = vars.image_file_name;
    sampling_pattern = vars.sampling_pattern;
    param_sampling = vars.param_sampling;
    im = vars.im;
    N = vars.N;
    Ny = vars.Ny;
    Nx = vars.Nx;
    uw = vars.uw;
    vw = vars.vw;
    u = vars.u;
    v = vars.v;
    uvidx = vars.uvidx;
    input_snr_v = vars.input_snr_v;

    %% compute weights
    param_precond_ = param_precond;
    
    param_precond_.N = N; % number of pixels in the image
    param_precond_.Nox = ox*Nx; % number of pixels in the image
    param_precond_.Noy = oy*Ny; % number of pixels in the image
    [aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond_);

    %% set the blocks structure

    if regenerate_block_structure
        [u_, v_, ~, uvidx_, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);

        R = length(u_);

        y0 = cell(num_tests, 1);
        y = cell(num_tests, 1);

        u = u_;
        v = v_;
        uvidx = uvidx_;

        u_ = [];
        v_ = [];
        uvidx_ = [];

        if use_real_visibilities
            for k = 1:num_tests
                y_ = cell(R, 1);
                for q = 1:R
                    y_{q} = yf{k}(uvidx{q});
                end
                y{k} = y_;
                y_ = [];
            end
        end

        if ~use_real_visibilities
            for k = 1:num_tests
                y0_ = cell(R, 1);
                y_ = cell(R, 1);
                for q = 1:R
                    y_{q} = yf{k}(uvidx{q});
                    y0_{q} = y0f{k}(uvidx{q});
                end
                y{k} = y_;
                y0{k} = y0_;
                y_ = []
                y0_ = [];
            end
        end
    end

    [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);

    wparam = [];
    wparam.N = N; % number of pixels in the image
    wparam.Nox = ox*Nx; % number of pixels in the image
    wparam.Noy = oy*Ny; % number of pixels in the image


    [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);


    %% compute the operator norm
    fprintf('Computing operator norms ...\n');

    fprintf('Natural W ...\n');
    evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);

    fprintf('No natural W ...\n');
    Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
    Gwt_ = Gw_';
    evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
    Gw_ = [];
    Gwt_ = [];        


    fprintf('Preconditioning ...\n');
    evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);

    evl_blocks = zeros(R, 1);

    if compute_block_op_norm 
        % maximum eigenvalue of operator
        for q = 1:R
            No = size(W{1}, 1);
            Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
            Tw_(:, W{q}) = T{q};
            fprintf('\nComputing operator norm: block %i \n', q)
            evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
            Tw_ = [];
        end
    end


    %% L2 ball sizes
    epsilon = cell(num_tests, 1);
    epsilons = cell(num_tests, 1);
    epsilonT = cell(num_tests, 1);
    epsilonTs = cell(num_tests, 1);
    %% generate bounds
    for k = 1:num_tests
        if use_real_visibilities
            [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
                [], sigma_noise, l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
        end
        if ~use_real_visibilities
            [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
                input_snr, [], l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
        end
    end
    %% gridding
    if use_gridded_data
        [yT, T, Tw] = util_compute_gridding_data(y, G, Gw);
    else
        T = G;
        Tw = Gw;
        yT = y;
    end

else
    error('Unable to find data files');
end
