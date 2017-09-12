if ~exist('results_prefix', 'var')
    results_prefix = '';
end


if save_data_on_disk
    % compute some data to save
    
    file_start = 0;
    while exist(sprintf('%s%s_%s_%s_%s_alg_data.mat', save_path, int2str(save_dataset_number), int2str(save_dataset_subnumber), int2str(file_start), results_prefix), 'file') || ... 
          exist(sprintf('%s%s_%s_%s_%s_alg_config.mat', save_path, int2str(save_dataset_number), int2str(save_dataset_subnumber), int2str(file_start), results_prefix), 'file')
        file_start = file_start + 1;
    end
    save(sprintf('%s%s_%s_%s_%s_alg_data.mat', save_path, int2str(save_dataset_number), int2str(save_dataset_subnumber), int2str(file_start), results_prefix), ...
        '-v7.3', 'result_st');
    try
        add_value = param_structure_name;
    catch
        add_value = '';
    end
    other = [];
    if ~isempty(who(sprintf('%s_*_vec', add_value)))
        other = [other sprintf('%s_*_vec', add_value)]; 
    end
    
    if ~use_real_visibilities
        save(sprintf('%s%s_%s_%s_%s_alg_config.mat', save_path, int2str(save_dataset_number), int2str(save_dataset_subnumber), int2str(file_start), results_prefix), ...
            'epsilonT', 'epsilonTs', 'epsilon', 'epsilons', 'num_tests', 'use_same_stop_criterion', sprintf('%s', add_value), ...
            'evl', 'evl_precond', 'evl_blocks', 'param_precond', 'aW', 'aWw', 'u', 'v', 'wlt_basis', 'nlevel', ...
            'param_block_structure', 'param_l2_ball', 'ox', 'oy', 'Kx', 'Ky', 'l2_ball_definition', 'stopping_criterion', ...
            'y0', 'y', 'uvidx');
    end
    if use_real_visibilities
        save(sprintf('%s%s_%s_%s_%s_alg_config.mat', save_path, int2str(save_dataset_number), int2str(save_dataset_subnumber), int2str(file_start), results_prefix), ...
            'epsilonT', 'epsilonTs', 'epsilon', 'epsilons', 'num_tests', 'use_same_stop_criterion', sprintf('%s', add_value), ...
            'evl', 'evl_precond', 'evl_blocks', 'param_precond', 'aW', 'aWw', 'u', 'v', 'wlt_basis', 'nlevel', ...
            'param_block_structure', 'param_l2_ball', 'ox', 'oy', 'Kx', 'Ky', 'l2_ball_definition', 'stopping_criterion', ...
            'y', 'uvidx');
    end
end