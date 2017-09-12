function util_create_pool(NumWorkers, cores)
% opens a new paralel pool with NumWorkers workers
% if a similar pool already exists it keeps it

global create_pool_name;
if isempty(create_pool_name)
    create_pool_name = 'local';
end

if ~isempty(getenv('SLURM_JOB_ID'))
    if ~isempty(gcp('nocreate'))
        location = strcat(strcat(getenv('SLURM_SUBMIT_DIR'), '/'), getenv('SLURM_JOB_ID'));
        try system(sprintf('mkdir %s', location)); end
        pc = parcluster(create_pool_name);
        if pc.NumWorkers ~= str2double(getenv('SLURM_CPUS_ON_NODE'))
            delete(gcp('nocreate'));
            pc = parcluster(create_pool_name);
            pc.JobStorageLocation = location;

            pc.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
            parpool(pc, str2double(getenv('SLURM_CPUS_ON_NODE')));

            fprintf('Created parallel pool ... \n');
        else
            fprintf('Using existing parallel pool ... \n');
        end
    else
        location = strcat(strcat(getenv('SLURM_SUBMIT_DIR'), '/'), getenv('SLURM_JOB_ID'));
        try system(sprintf('mkdir %s', location)); end
        pc = parcluster(create_pool_name);
        pc.JobStorageLocation = location;

        pc.NumWorkers = str2double(getenv('SLURM_CPUS_ON_NODE'));
        parpool(pc, str2double(getenv('SLURM_CPUS_ON_NODE')));

        fprintf('Created parallel pool ... \n');
    end
else
    if ~isempty(gcp('nocreate'))
        pc = parcluster(create_pool_name);
        if pc.NumWorkers ~= NumWorkers
            delete(gcp('nocreate'));
            pc = parcluster(create_pool_name);
            pc.NumWorkers = NumWorkers;
            parpool(pc, NumWorkers);
            fprintf('Created parallel pool ... \n');
        else
            fprintf('Using existing parallel pool ... \n');
        end
    else
        pc = parcluster(create_pool_name);
        pc.NumWorkers = NumWorkers;
        parpool(pc, NumWorkers);
        fprintf('Created parallel pool ... \n');
    end
end

fprintf('Setting DWT mode ... \n');
spmd
    dwtmode('per');
end

end

