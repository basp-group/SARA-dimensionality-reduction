% script that loads or generates the input data for the test


% set random nr gen
rng shuffle;


if gen_data == 0
    fprintf('Loading data from disk ... \n\n');
    if exist(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)), 'file') || ... 
          exist(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)), 'file')
      
        load(sprintf('%s%s_input_data.mat', save_path, int2str(save_dataset_number)));
        load(sprintf('%s%s_input_data_config.mat', save_path, int2str(save_dataset_number)));

        
        %% compute weights
        param_precond.N = N; % number of pixels in the image
        param_precond.Nox = ox*Nx; % number of pixels in the image
        param_precond.Noy = oy*Ny; % number of pixels in the image
        [aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);
        
        %% set the blocks structure

        if regenerate_block_structure
            [u_, v_, ~, uvidx_, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
        
            R = length(u_);
            
            y0 = cell(num_tests, 1);
            y = cell(num_tests, 1);

            u = u_;
            v = v_;
            uvidx = uvidx_;
            
            clear u_ v_ uvidx_;
            
            if use_real_visibilities
                for k = 1:num_tests
                    y_ = cell(R, 1);
                    for q = 1:R
                        y_{q} = yf{k}(uvidx{q});
                    end
                    y{k} = y_;
                    clear y_;
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
                    clear y_ y0_;
                end
            end
        else
            R = length(u);
            nW = cell(R, 1);
            aW = cell(R, 1);
            for q = 1:R
                nW{q} = nWw(uvidx{q});
                aW{q} = aWw(uvidx{q});
            end
        end
        %%
        [A, At, G, W, Gw, As, Ats, S] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, param_nufft);
        
        wparam.N = N; % number of pixels in the image
        wparam.Nox = ox*Nx; % number of pixels in the image
        wparam.Noy = oy*Ny; % number of pixels in the image
        
        
        [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        
        
        %% compute the operator norm
        fprintf('Computing operator norms ...\n');
        
        if compute_evl
            fprintf('Natural W ...\n');
            evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping natural W ...\n');
            evl = 0;
        end

        if compute_evl_no_natw
            fprintf('No natural W ...\n');
            Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
            Gwt_ = Gw_';
            evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
            clear Gw_ Gwt_;
        else
            fprintf('Skipping no natural W ...\n');
            evl_no_natw = 0;
        end
        
        if compute_evl_precond
            fprintf('Preconditioning ...\n');
            evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping preconditioning ...\n');
            evl_precond = 0;
        end
        
        
        evl_blocks = zeros(R, 1);

        if compute_block_op_norm 
            % maximum eigenvalue of operator
            for q = 1:R
                No = size(W{1}, 1);
                Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
                Tw_(:, W{q}) = T{q};
                fprintf('\nComputing operator norm: block %i \n', q)
                evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
                clear Tw_;
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
end

if gen_data == 1
    
    if use_real_visibilities
        fprintf('Using real data ... \n\n');
        
        Nx = param_real_data.image_size_Nx;
        Ny = param_real_data.image_size_Ny;
        N = Nx * Ny;
        
        %% generate the sampling pattern
        param_sampling.N = N; % number of pixels in the image
        param_sampling.Nox = ox*Nx; % number of pixels in the image
        param_sampling.Noy = oy*Ny; % number of pixels in the image

        %% load noisy input data
        yf = cell(num_tests, 1);
        y0f = cell(num_tests, 1);
        
%         [yf{num_tests}, uw, vw, sigma_noise, nWw] = util_load_real_vis_data(visibility_file_name, param_real_data);
        [yf{num_tests}, uw, vw, sigma_noise, nWw] = util_load_hyperspectral_averaged_data(visibility_file_name, param_real_data);
        if use_symmetric_fourier_sampling
            uw = [uw; -uw];
            vw = [vw; -vw];
            yf{num_tests} = [yf{num_tests}(:); conj(yf{num_tests}(:))];
            nWw = [nWw; nWw];
        end
        
        %% compute weights
        param_precond.N = N; % number of pixels in the image
        param_precond.Nox = ox*Nx; % number of pixels in the image
        param_precond.Noy = oy*Ny; % number of pixels in the image
        [aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);
        
        %% set the blocks structure
        [u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);

        
        
        %% measurement operator initialization 

        fprintf('Initializing the NUFFT operator\n\n');
        tstart = tic;

        [A, At, G, W, Gw, As, Ats, S] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, param_nufft);
        
        tend = toc(tstart);

        fprintf('Initialization runtime: %ds\n\n', ceil(tend));
        R = length(v);

        y0 = cell(num_tests, 1);
        y = cell(num_tests, 1);
        
        epsilon = cell(num_tests, 1);
        epsilons = cell(num_tests, 1);
        epsilonT = cell(num_tests, 1);
        epsilonTs = cell(num_tests, 1);

        for k = 1:num_tests
            y_ = cell(R, 1);
            for q = 1:R
                y_{q} = yf{k}(uvidx{q});
            end
            y{k} = y_;
            clear y_;
            [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
                [], sigma_noise, l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
        end

        %% gridding
        if use_gridded_data
            [yT, T, Tw] = util_compute_gridding_data(y, G, Gw);
        else
            T = G;
            Tw = Gw;
            yT = y;
        end

        %% sparsity operator definition

        [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);


        %% compute the operator norm
        fprintf('Computing operator norms ...\n');
        
        if compute_evl
            fprintf('Natural W ...\n');
            evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping natural W ...\n');
            evl = 0;
        end

        if compute_evl_no_natw
            fprintf('No natural W ...\n');
            Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
            Gwt_ = Gw_';
            evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
            clear Gw_ Gwt_;
        else
            fprintf('Skipping no natural W ...\n');
            evl_no_natw = 0;
        end
        
        if compute_evl_precond
            fprintf('Preconditioning ...\n');
            evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping preconditioning ...\n');
            evl_precond = 0;
        end
        
        
        evl_blocks = zeros(R, 1);

        if compute_block_op_norm 
            % maximum eigenvalue of operator
            for q = 1:R
                No = size(W{1}, 1);
                Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
                Tw_(:, W{q}) = T{q};
                fprintf('\nComputing operator norm: block %i \n', q)
                evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
                clear Tw_;
            end
        end

        %% save data
        if save_data_on_disk == 1
            fprintf('Saving new data ... \n');

            if save_data_on_disk
                file_start = save_dataset_number;
                while exist(sprintf('%s%s_input_data.mat', save_path, int2str(file_start)), 'file') || ... 
                      exist(sprintf('$s%s_input_data_config.mat', save_path, int2str(file_start)), 'file')
                    file_start = file_start + 1;
                end
                if file_start ~= save_dataset_number
                    fprintf('WARNING: Saving new data in file %d instead of %d \n\n', file_start, save_dataset_number);
                end

                save(sprintf('%s%s_input_data', save_path, int2str(file_start)), '-v7.3', 'y', ... % 'G', 'Gw', 'W'
                    'yf', 'nWw'); % do not save G, it is large and is faster to compute than saving 
                save(sprintf('%s%s_input_data_config', save_path, int2str(file_start)), 'N', 'Ny', ...
                    'Nx', 'uw', 'vw', 'u', 'v', 'uvidx', 'param_sampling', 'sigma_noise', 'visibility_file_name', 'param_real_data', 'num_tests', 'use_real_visibilities');
            end
        end
    end
    
    
    if use_simulated_data
        fprintf('Generating new data ... \n\n');
        %% image loading
        [im, N, Ny, Nx] = util_read_image(image_file_name);

        %% generate the sampling pattern
        param_sampling.N = N; % number of pixels in the image
        param_sampling.Nox = ox*Nx; % number of pixels in the image
        param_sampling.Noy = oy*Ny; % number of pixels in the image
        
        [uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);

        
        if use_symmetric_fourier_sampling
            uw = [uw; -uw];
            vw = [vw; -vw];
        end
        
        %% compute weights
        param_precond.N = N; % number of pixels in the image
        param_precond.Nox = ox*Nx; % number of pixels in the image
        param_precond.Noy = oy*Ny; % number of pixels in the image
        [aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);
        
        %% set the blocks structure
        nWw = ones(length(uw), 1);
        [u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);

        %% measurement operator initialization 

        fprintf('Initializing the NUFFT operator\n\n');
        tstart = tic;
        [A, At, G, W, Gw, As, Ats, S] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, param_nufft);
        
        tend = toc(tstart);

        fprintf('Initialization runtime: %ds\n\n', ceil(tend));
        R = length(v);

        y0 = cell(num_tests, 1);
        y = cell(num_tests, 1);
        aY = cell(num_tests, 1);
        epsilon = cell(num_tests, 1);
        epsilons = cell(num_tests, 1);
        epsilonT = cell(num_tests, 1);
        epsilonTs = cell(num_tests, 1);
        y0f = cell(num_tests, 1);
        yf = cell(num_tests, 1);
        input_snr_v = cell(num_tests, 1);
        
        %% generate noisy input data
        for k = 1:num_tests
            [y0{k}, y{k}, y0f{k}, yf{k}, aY{k}, input_snr_v{k}, im_C] = util_gen_input_data(im, G, W, A, input_snr, ...
                use_different_per_block_input_snr, per_block_input_snr_delta, uvidx, ...
                param_image_var);

            if use_symmetric_fourier_sampling
                y0f{k} = [y0f{k}(uvidx{k}(1:end/2)); conj(y0f{k}(uvidx{k}(1:end/2)))];
                yf{k} = [yf{k}(uvidx{k}(1:end/2)); conj(yf{k}(uvidx{k}(1:end/2)))];
                for j = 1:R
                    y{k}{j} = [y{k}{j}(uvidx{k}(1:end/2)); conj(y{k}{j}(uvidx{k}(1:end/2)))];
                    y0{k}{j} = [y0{k}{j}(uvidx{k}(1:end/2)); conj(y0{k}{j}(uvidx{k}(1:end/2)))];
                    aY{k}{j} = [aY{k}{j}(uvidx{k}(1:end/2)); conj(aY{k}{j}(uvidx{k}(1:end/2)))];
                end
            end
            
            
            [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
                input_snr, [], l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
        end

        %% gridding
        if use_gridded_data
            [yT, T, Tw] = util_compute_gridding_data(y, G, Gw);
        else
            T = G;
            Tw = Gw;
            yT = y;
        end

        %% sparsity operator definition

        [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);


        %% compute the operator norm
        fprintf('Computing operator norms ...\n');
        
        if compute_evl
            fprintf('Natural W ...\n');
            evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping natural W ...\n');
            evl = 0;
        end

        if compute_evl_no_natw
            fprintf('No natural W ...\n');
            Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
            Gwt_ = Gw_';
            evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
            clear Gw_ Gwt_;
        else
            fprintf('Skipping no natural W ...\n');
            evl_no_natw = 0;
        end
        
        if compute_evl_precond
            fprintf('Preconditioning ...\n');
            evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
        else
            fprintf('Skipping preconditioning ...\n');
            evl_precond = 0;
        end
        
        
        evl_blocks = zeros(R, 1);

        if compute_block_op_norm 
            % maximum eigenvalue of operator
            for q = 1:R
                No = size(W{1}, 1);
                Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
                Tw_(:, W{q}) = T{q};
                fprintf('\nComputing operator norm: block %i \n', q)
                evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
                clear Tw_;
            end
        end

        %% save data
        if save_data_on_disk == 1
            fprintf('Saving new data ... \n');

            if save_data_on_disk
                file_start = save_dataset_number;
                while exist(sprintf('%s%s_input_data.mat', save_path, int2str(file_start)), 'file') || ... 
                      exist(sprintf('$s%s_input_data_config.mat', save_path, int2str(file_start)), 'file')
                    file_start = file_start + 1;
                end
                if file_start ~= save_dataset_number;
                    fprintf('WARNING: Saving new data in file %d instead of %d \n\n', file_start, save_dataset_number);
                end

                save(sprintf('%s%s_input_data', save_path, int2str(file_start)), '-v7.3', ... % 'G', 'Gw', 'W'
                    'y0f', 'yf', 'y', 'y0', 'nWw'); % do not save G, it is large and is faster to compute than saving 
                save(sprintf('%s%s_input_data_config', save_path, int2str(file_start)), 'N', 'Ny', ...
                    'Nx', 'uw', 'vw', 'u', 'v', 'uvidx', 'im', 'im_C', 'sampling_pattern', 'param_sampling', 'input_snr', 'input_snr_v', 'image_file_name', 'num_tests', 'use_real_visibilities');

                if strcmp(sampling_pattern, 'file')
                    for k = 1:num_tests
                        vis = yf{k};
                        save(sprintf('%s%s_vis_data_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                        clear vis;
                        
                        vis = yf{k} - y0f{k};
                        save(sprintf('%s%s_vis_data_noise_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                        clear vis;
                    end
                end
            end
        end
    end
end

if gen_data == 2
    fprintf('Using data from workspace ... \n\n');
end

if free_memory
    % free memory of the whose Gw or the split G are not needed
    try 
        if ~run_admm_bpcon && ~run_sdmm_bpcon && ~run_fb_nnls && ~run_krylov_nnls
            clear Gw;
            clear Tw;
        end

        if ~run_pdfb_bpcon_par_sim_rescaled_precond_wave_par && ~run_pdfb_bpcon_par_sim && ~run_pdfb_bpcon_par_sim_rescaled && ~run_pdfb_bpcon_par_sim_rescaled_precond &&... 
            ~run_pdfb_bpcon_par_sim_rand_rescaled && ~run_pdfb_bpcon_par_sim_block_rand_rescaled && ...
            ~run_pdfb_bpcon_dist && ~run_pdfb_bpcon_dist_rescaled && ~run_admm_bpconpar && ~run_sdmm_bpconpar && ...
            ~run_pdfb_bpcon_par_sim_rescaled_rec_async && ...
            ~run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p && ...
            ~run_admm_bpconpar_wavepar && ~run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps && ...
            ~run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps
            clear G;
            clear T;
            clear W;
        end
    end
end
