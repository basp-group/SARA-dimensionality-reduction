util_create_pool(num_workers);

if run_pdfb_bpcon_par_sim
    
    result_st = [];

    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, L2_vp{i}] ...
            = pdfb_bpcon_par_sim(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

if run_pdfb_bpcon_par_sim_rescaled
    
    result_st = [];
     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, ...
            L2_vp{i}, delta_v{i}, sol_v{i}, snr_v{i}, ~, ~, sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, W, Psi, Psit, Psiw, Psitw, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.sol_reweight_v = sol_reweight_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end



if run_pdfb_bpcon_par_sim_rescaled_precond
    
    result_st = [];
     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, ...
            L2_vp{i}, delta_v{i}, sol_v{i}, snr_v{i}, no_sub_itr_v{i}, ~, ~, sol_best_bound_v{i}, sol_reweight_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_precond(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_precond runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.sol_best_bound_v = sol_best_bound_v;
    result_st.sol_reweight_v = sol_reweight_v;
    result_st.snr_v = snr_v;
    result_st.no_sub_itr_v = no_sub_itr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond';
    param_structure_name = 'param_pdfb_precond';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end



if run_pdfb_bpcon_par_sim_rescaled_rec_async
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rescaled_rec_async\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, ...
            L2_vp{i}, delta_v{i}, sol_v{i}, snr_v{i}] ...
            = pdfb_bpcon_par_sim_rescaled_rec_async(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_oas);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rescaled_rec_async runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_rescaled_rec_async';
    param_structure_name = 'param_pdfb_oas';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

    
if run_pdfb_bpcon_par_sim_rand_rescaled
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rand_rescaled\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, ...
            L2_vp{i}, delta_v{i}, sol_v{i}, snr_v{i}] ...
            = pdfb_bpcon_par_sim_rand_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, W, Psi, Psit, param_pdfb_prob);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_rand_rescaled';
    param_structure_name = 'param_pdfb_prob';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

    
if run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, L2_vp{i}, ...
            delta_v{i}, sol_v{i}, snr_v{i}] = pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_prob_ad);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p';
    param_structure_name = 'param_pdfb_prob_ad';

    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

if run_pdfb_bpcon_par_sim_block_rand_rescaled
    result_st = [];

    
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_par_sim_block_rand_rescaled\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, L2_vp{i}] ...
            = pdfb_bpcon_par_sim_block_rand_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb_prob_bc);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_par_sim_block_rand_rescaled runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_par_sim_block_rand_rescaled';
    param_structure_name = 'param_pdfb_prob_bc';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

if run_pdfb_bpcon_dist
    result_st = [];

    
    parfor i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_dist\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, L2_vp{i}] ...
            = pdfb_bpcon_dist(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_dist runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_dist';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

if run_pdfb_bpcon_dist_rescaled
    result_st = [];

    
    parfor i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running pdfb_bpcon_dist_rescaled\n');
        [sol{i}, L1_v{i}, L1_vp{i}, L2_v{i}, L2_vp{i}] ...
            = pdfb_bpcon_dist_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, A, At, T, W, Psi, Psit, param_pdfb);
        tend = toc(tstart_a);
        fprintf(' pdfb_bpcon_dist_rescaled runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L1_vp = L1_vp;
    result_st.L2_v = L2_v;
    result_st.L2_vp = L2_vp;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'pdfb_bpcon_dist_rescaled';
    param_structure_name = 'param_pdfb';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


if run_admm_bpcon
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpcon\n');
        [sol{i}, ~, L1_v{i}, L2_v{i}] = admm_bpcon(cell2mat(yT{i}), epsilon{i}, epsilons{i}, @(x) Tw * A(x), @(x) At(Tw' * x), Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpcon runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L2_v = L2_v;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'admm_bpcon';
    param_structure_name = 'param_admm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


if run_admm_bpconpar
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running admm_bpconpar\n');
        [sol{i}, ~, L1_v{i}, L2_v{i}, ...
            delta_v{i}, sol_v{i}, snr_v{i}] ...
            = admm_bpconpar(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T, W, Psiw, Psitw, param_admm);
        tend = toc(tstart_a);
        fprintf(' admm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L2_v = L2_v;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'admm_bpconpar';
    param_structure_name = 'param_admm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


if run_sdmm_bpcon
    result_st = [];

     
    parfor i = 1:num_tests
        % wavelet mode is a global variable which does not get transfered
        % to the workes; we need to set it manually for each worker
        dwtmode('per');

        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running sdmm_bpcon\n');
        [sol{i}, L1_v{i}, L2_v{i}] = sdmm_bpcon(cell2mat(yT{i})/sqrt(evl), epsilon{i}/sqrt(evl), epsilons{i}/sqrt(evl), @(x) Tw * A(x)/sqrt(evl), @(x) At(Tw' * x)/sqrt(evl), Psiw, Psitw, param_sdmm);
        L1_v{i} = L1_v{i} * sqrt(evl);
        L2_v{i} = L2_v{i} * sqrt(evl);
        sol{i} = real(sol{i});
        sol{i}(sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' sdmm_bpcon runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L2_v = L2_v;
    result_st.time = time;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'sdmm_bpcon';
    param_structure_name = 'param_sdmm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


if run_sdmm_bpconpar
    result_st = [];

     
    parfor i = 1:num_tests
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
        [sol{i}, L1_v{i}, L2_v{i}, ...
            delta_v{i}, sol_v{i}, snr_v{i}] ...
            = sdmm_bpconpar(yT_, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), @(x) A(x), @(x) At(x), T, W, Psiw, Psitw, param_sdmm);
        sol{i} = real(sol{i});
        sol{i}(sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' sdmm_bpconpar runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L1_v = L1_v;
    result_st.L2_v = L2_v;
    result_st.time = time;
    result_st.delta_v = delta_v;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'sdmm_bpconpar';
    param_structure_name = 'param_sdmm';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end

if run_fb_nnls
    result_st = [];
    
    
    parfor i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running fb_nnls\n');
        [sol{i}, L2_v{i}, sol_v{i}, snr_v{i}] = fb_nnls(cell2mat(yT{i}), @(x) Tw * A(x), @(x) At(Tw' * x), param_nnls);
        sol{i} = real(sol{i});
        sol{i}(sol{i}<0) = 0;
        tend = toc(tstart_a);
        fprintf(' fb_nnls runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L2_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.L2_v = L2_v;
    result_st.time = time;
    result_st.sol_v = sol_v;
    result_st.snr_v = snr_v;
        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'nnls';
    param_structure_name = 'param_nnls';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


if run_krylov_nnls
    result_st = [];

    parfor i = 1:num_tests
        fprintf('Test run %i:\n', i);
    
        tstart_a = tic;
        fprintf(' Running krylov_nnls\n');
        [sol__] = fcgls_rp_nn_ReSt(complex2vec(cell2mat(yT{i})), @(x) complex2vec(Tw * A(reshape((x), Ny, Nx))), @(x) real(reshape(At(Tw' * vec2complex(x)), Ny*Nx, 1)), param_krylov_nnls);
        sol{i} = reshape(sol__(:, end), Ny, Nx);
        
        tend = toc(tstart_a);
        fprintf(' krylov_nnls runtime: %ds\n\n', ceil(tend));
        
        time{i} = tend;
        error = im - sol{i};
        snr_end{i} = 20 * log10(norm(im(:))/norm(error(:)));
        no_itr{i} = length(L1_v{i});

        wcoef = [];
        for q = 1:length(Psit)
            wcoef = [wcoef; Psit{q}(sol{i})];
        end
        sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
    end
    
    result_st.sol = sol;
    result_st.time = time;

        
    result_st.snr = snr_end;
    result_st.sparsity = sparsity;
    result_st.no_itr = no_itr;
    
    results_prefix = 'krylov_nnls';
    param_structure_name = 'param_krylov_nnls';
    
    % results
    script_gen_figures;

    % save data
    script_save_result_data;
end


