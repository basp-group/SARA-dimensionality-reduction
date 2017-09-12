%% PD with simulated parallelism

if run_pdfb_bpcon_par_sim
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, result_st.L2_vp{i}] ...
            = pdfb_bpcon_par_sim(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism

if run_pdfb_bpcon_par_sim_rescaled
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled\n');
%         [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
%             result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, ~, ~, result_st.sol_reweight_v{i}] ...
%             = pdfb_bpcon_par_sim_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, B, Bt, T, W, Psi, Psit, Psiw, Psitw, param_pdfb);
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, ~, ~, result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_adapt_eps(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, B, Bt, T, W, Psi, Psit, Psiw, Psitw, param_pdfb);

        
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescalled simulated PD with simulated parallelism and natural weigths 
%% in the projection instead of the operator

if run_pdfb_bpcon_par_sim_rescaled_natw
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

    fprintf('The natural weighting is included in the operator G \n\n Undoing natural weighting in the operator G\n');
    T_ = cell(size(T));
    for q = 1:length(T)
        T_{q} = spdiags(1./nW{q}, 0, length(nW{q}), length(nW{q})) * T{q};
    end
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_natw\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, ~, ~, result_st.sol_reweight_v{i}, result_st.no_sub_itr_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_natw(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T_, W, Psi, Psit, Psiw, Psitw, nW, param_pdfb_natw);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_natw runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    clear T_;
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_natw';
    param_structure_name = 'param_pdfb_natw';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
if run_pdfb_bpcon_par_sim_rescaled_precond
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.sol_best_bound_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
    result_st.xcorr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
            result_st.sol_best_bound_v{i}, result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond';
    param_structure_name = 'param_pdfb_precond';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% parallelised wavelets
if run_pdfb_bpcon_par_sim_rescaled_precond_wave_par
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.sol_best_bound_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
    result_st.xcorr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond_wave_par\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
            result_st.sol_best_bound_v{i}, result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond_wave_par(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond_wave_par runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond_wave_par';
    param_structure_name = 'param_pdfb_precond';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% parallelised wavelets and adaptive eps and gamma
if run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    result_st.v1 = cell(num_tests, 1);
    result_st.v2 = cell(num_tests, 1);

    if free_memory
        G = [];
        uw = [];
        vw = [];
        v = [];
        u = [];
        nWw = [];
        nW = [];
        im_C = [];
        im = [];
        aWw = [];
    end
    
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ...
            result_st.v1{i}, result_st.v2{i}, ...
            result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond_eps);
        tend = toc(tstart_a);
        fprintf(' run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps';
    param_structure_name = 'param_pdfb_precond_eps';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% adaptive eps and gamma
if run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.sol_reweight_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    result_st.v1 = cell(num_tests, 1);
    result_st.v2 = cell(num_tests, 1);

    if free_memory
        G = [];
        uw = [];
        vw = [];
        v = [];
        u = [];
        nWw = [];
        nW = [];
        im_C = [];
        im = [];
        aWw = [];
    end
    
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ...
            result_st.v1{i}, result_st.v2{i}, ...
            result_st.sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond_var_block_eps(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond_eps);
        tend = toc(tstart_a);
        fprintf(' run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps';
    param_structure_name = 'param_pdfb_precond_eps';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end



%% rescaled PD with simulated parallelism and preconditioning
%% parallelised wavelets
%% golden search for best epsilon
if run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_gs
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond_wave_par_gs\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond_wave_par_golden_search(yT{i}, sigma_noise, As, Ats, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond_gs);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond_wave_par_gs runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond_wave_par_gs';
    param_structure_name = 'param_pdfb_precond_gs';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% GPU fft
if run_pdfb_bpcon_par_sim_rescaled_gpu
    
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);

     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_gpu\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_gpu(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_gpu runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
if run_pdfb_bpcon_par_sim_rescaled_rec_async
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_rec_async\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_rec_async(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_oas);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_rec_async runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_rec_async';
    param_structure_name = 'param_pdfb_oas';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end
    
%% rescaled PD with simulated parallelism and preconditioning
%% randomised
if run_pdfb_bpcon_par_sim_rand_rescaled
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rand_rescaled\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
            result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = pdfb_bpcon_par_sim_rand_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, W, Psi, Psit, param_pdfb_prob);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'pdfb_bpcon_par_sim_rand_rescaled';
    param_structure_name = 'param_pdfb_prob';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% randomised with non uniform probabilities
if run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, result_st.L2_vp{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] = pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_prob_ad);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p';
    param_structure_name = 'param_pdfb_prob_ad';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% rescaled PD with simulated parallelism and preconditioning
%% randomised with non uniform probabilities
if run_pdfb_bpcon_par_sim_block_rand_rescaled
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_block_rand_rescaled\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, result_st.L2_vp{i}] ...
            = pdfb_bpcon_par_sim_block_rand_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_prob_bc);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_block_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    results_prefix = 'pdfb_bpcon_par_sim_block_rand_rescaled';
    param_structure_name = 'param_pdfb_prob_bc';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% basic PD 
if run_pdfb_bpcon_dist
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    
    
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    
    for i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_dist\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, result_st.L2_vp{i}] ...
            = pdfb_bpcon_dist(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_dist runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'pdfb_bpcon_dist';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% basic PD with distributed data fidelity
if run_pdfb_bpcon_dist_rescaled
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    
    
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    
    for i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_dist_rescaled\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, result_st.L2_vp{i}] ...
            = pdfb_bpcon_dist_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_dist_rescaled runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'pdfb_bpcon_dist_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% ADMM

if run_admm_bpcon
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);


    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpcon\n');
        [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}] = admm_bpcon(cell2mat(yT{i}), epsilon{i}, epsilons{i}, @(x) Tw * A(x), @(x) At(Tw' * x), Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpcon runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'admm_bpcon';
    param_structure_name = 'param_admm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% ADMM

if run_admm_bpconpar
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpconpar\n');
        [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = admm_bpconpar(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T, W, Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'admm_bpconpar';
    param_structure_name = 'param_admm';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% ADMM
if run_admm_bpconpar_natw
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    result_st.no_sub_itr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    
    fprintf('The natural weighting is included in the operator G \n\n Undoing natural weighting in the operator G\n');
    T_ = cell(size(T));
    for q = 1:length(T)
        T_{q} = spdiags(1./nW{q}, 0, length(nW{q}), length(nW{q})) * T{q};
    end
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpconpar_natw\n');
        [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}] ...
            = admm_bpconpar_natw(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T_, W, Psiw, Psitw, nW, param_admm_natw);
        tend = toc(tstart_a);
        fprintf(' admm_bpconpar_natw runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    clear T_;
    
    results_prefix = 'admm_bpconpar_natw';
    param_structure_name = 'param_admm_natw';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% ADMM with parallel wavelets
if run_admm_bpconpar_wavepar
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpconpar_wavepar\n');
        [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = admm_bpconpar_wavepar(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T, W, Psi, Psit, Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'admm_bpconpar_wavepar';
    param_structure_name = 'param_admm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% SDMM
if run_sdmm_bpcon
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);


    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running sdmm_bpcon\n');
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L2_v{i}] = sdmm_bpcon(cell2mat(yT{i})/sqrt(evl), epsilon{i}/sqrt(evl), epsilons{i}/sqrt(evl), @(x) Tw * A(x)/sqrt(evl), @(x) At(Tw' * x)/sqrt(evl), Psiw, Psitw, param_sdmm);
        result_st.L1_v{i} = result_st.L1_v{i} * sqrt(evl);
        result_st.L2_v{i} = result_st.L2_v{i} * sqrt(evl);
        result_st.sol{i} = real(sol{i});
        result_st.sol{i}(result_st.sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' sdmm_bpcon runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'sdmm_bpcon';
    param_structure_name = 'param_sdmm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% SDMM
if run_sdmm_bpconpar
    result_st = [];
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
    result_st.delta_v = cell(num_tests, 1);
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
     
    for i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running sdmm_bpconpar\n');
        yT_ = cell(length(yT{i}), 1);
        for o=1:length(yT{i})
            yT_{o} = yT{i}{o};
        end
        [result_st.sol{i}, result_st.L1_v{i}, result_st.L2_v{i}, ...
            result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
            = sdmm_bpconpar(yT_, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), @(x) A(x), @(x) At(x), T, W, Psiw, Psitw, param_sdmm);
        result_st.sol{i} = real(result_st.sol{i});
        result_st.sol{i}(result_st.sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' sdmm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'sdmm_bpconpar';
    param_structure_name = 'param_sdmm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% NNLS
if run_fb_nnls
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.L1_v = cell(num_tests, 1);
    result_st.L1_vp = cell(num_tests, 1);
    result_st.L2_v = cell(num_tests, 1);
    result_st.L2_vp = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);
     
    result_st.sol_v = cell(num_tests, 1);
    result_st.snr_v = cell(num_tests, 1);
    
    for i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running fb_nnls\n');
        [result_st.sol{i}, result_st.L2_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] = fb_nnls(cell2mat(yT{i}), B, Bt, param_nnls);
        result_st.sol{i} = real(result_st.sol{i});
        result_st.sol{i}(result_st.sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' fb_nnls runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L2_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'nnls';
    param_structure_name = 'param_nnls';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

%% krylov tests
if run_krylov_nnls
    result_st = [];
    
    result_st.sol = cell(num_tests, 1);
    result_st.time = cell(num_tests, 1);

        
    result_st.snr = cell(num_tests, 1);
    result_st.sparsity = cell(num_tests, 1);
    result_st.no_itr = cell(num_tests, 1);
    
    for i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running krylov_nnls\n');
        [sol__] = fcgls_rp_nn_ReSt(complex2vec(cell2mat(yT{i})), @(x) complex2vec(Tw * A(reshape((x), Ny, Nx))), @(x) real(reshape(At(Tw' * vec2complex(x)), Ny*Nx, 1)), param_krylov_nnls);
        result_st.sol{i} = reshape(sol__(:, end), Ny, Nx);
        
        tend = toc(tstart_a);
        fprintf(' krylov_nnls runtime: %ds\n\n', ceil(tend));
        
        result_st.time{i} = tend;
        if ~use_real_visibilities
            error = im - result_st.sol{i};
            result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
        end
        result_st.no_itr{i} = length(result_st.L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(result_st.sol{i})];
        end
        result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    results_prefix = 'krylov_nnls';
    param_structure_name = 'param_krylov_nnls';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


