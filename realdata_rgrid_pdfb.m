addpath data/
addpath data/vis
addpath lib/

try
    run('src/irt/setup.m');
catch ME
    error('NUFFT library not found in location src/irt');
end

%% run parameters
% 0 - loads new data from file based on the dataset number supplied
% 1 - generates new data
% 2 - uses the data in matlab's workspace
gen_data = 1;
gen_figures = 1;
gen_only_average_figures = 0;
free_memory = 0;

save_dataset_number = 6; % number of the dataset to write files to
save_dataset_subnumber = 0; % number of the dataset to write files to

save_data_on_disk = 0; % flag
save_eps_files = 0; % flag
save_path = 'results/rgrid/';

num_tests = 1;
num_workers = 1; % number of tests to run in parallel; should be less than 
                 % the number of cores and is limited by the system memory for the variables

run_pdfb_bpcon_par_sim_rescaled = 1; % flag
run_pdfb_bpcon_par_sim_rescaled_natw = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_precond = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_precond_wave_par = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_gs = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_precond_wave_par_var_block_eps = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_precond_var_block_eps = 0;
run_pdfb_bpcon_par_sim_rand_rescaled = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_gpu = 0; % flag
run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p = 0; % flag

run_admm_bpconpar = 0; %flag
run_admm_bpconpar_natw = 0; %flag
run_admm_bpconpar_wavepar = 0; %flag
run_sdmm_bpconpar = 0; %flag

% work in progress algos
run_pdfb_bpcon_par_sim = 0; % flag
run_pdfb_bpcon_par_sim_block_rand_rescaled = 0; % flag
run_pdfb_bpcon_dist = 0; % flag
run_pdfb_bpcon_dist_rescaled = 0; % flag
run_pdfb_bpcon_par_sim_rescaled_rec_async = 0; % flag
run_krylov_nnls = 0; %flag

% old algos
run_admm_bpcon = 0; % flag
run_sdmm_bpcon = 0; % flag

run_fb_nnls = 0; % flag

%% real data generation
use_real_visibilities = 1;
visibility_file_name = 'data/vis/WEIGHTED-CYGA-C-6680-64CH';

param_real_data.image_size_Nx = 256;
param_real_data.image_size_Ny = 256;
param_real_data.pixel_size = 2.5;
param_real_data.use_shift = 0;
param_real_data.use_undersamplig = 0;

if use_real_visibilities % force only one test
    num_tests = 1;
    num_workers = 1;
end

%% simulated data generation
use_simulated_data = 0;

%% various config parameters
verbosity = 1;

nlevel = 4; % wavelet level
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

use_gridded_data = 0; % flag setting for generating gridded data

% evl params

compute_evl = 0;
compute_evl_no_natw = 0;
compute_evl_precond = 0;
compute_block_op_norm = 0; % flag to compute the operator norm for each block

use_symmetric_fourier_sampling = 0;

%% definition for the stopping criterion
% options: 
% l2_ball_definition -> 'sigma', 'chi-percentile', 'value'
% stopping_criterion -> 'sigma', 'chi-percentile', 'l2-ball-percentage', 'value'

l2_ball_definition = 'value';
stopping_criterion = 'l2-ball-percentage';

param_l2_ball.stop_eps_v = step_epsilon; % set epsilon value BEFORE running this script
param_l2_ball.val_eps_v = 1.0*param_l2_ball.stop_eps_v;

param_l2_ball.sigma_ball = 2;
param_l2_ball.sigma_stop = 2;

param_l2_ball.chi_percentile_ball = 0.99;
param_l2_ball.chi_percentile_stop = 0.999;

param_l2_ball.l2_ball_percentage_stop = 1.0001;

use_same_stop_criterion = 1; % forces the distributed criterion to be scaled
                             % such that same norm is imposed as in the nondistributed setup
 
%% sparsity prior
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used

%% nufft parameters

param_nufft.gen_fft_op_without_scale = 0;
param_nufft.use_fft_mask = 0;
param_nufft.use_fft_on_gpu = 0; % gpu FFT
param_nufft.use_nufft_blocks = 0;

%% block structure

regenerate_block_structure = 1;
param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;
param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;
param_block_structure.use_manual_frequency_partitioning = 0;
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 0;
param_block_structure.partition = [1000 2000 4000];

param_block_structure.use_equal_partitioning = 1;
param_block_structure.equal_partitioning_no = 1;

%% preconditioning

param_precond.gen_uniform_weight_matrix = 0; %set weighting type
param_precond.uniform_weight_sub_pixels = 1;

%% get input data
script_get_input_data;

%% For dimensionality reduction
% Flags

% Flag to pull up the values of elements of the holographic matrix
% This is to avoid having VERY small values which might later explode
% during computation of inverse or reciprocal.
thresholdholographic = 1; diagthresholdepsilon = 1e-10;

% Flag to set if we want to approximate D with an
% identity matrix. Reset the flag to use the normal D.
% (D is the diagonal aproximation of the covariance matrix)
identityapprox = 0;

% Percent of (ox*oy *image size) to keep after dimensionality reduction
klargestpercent = 100;

% Compute holographic matrix
h = Gw'*Gw;

% Create the new measurement operator
serialise = @(x) x(:); 
grid_fwd = @(x) h*A(x); % Phi^TPhi; input = [Ny, Nx] image; output = [Ny, Nx] matrix
grid_adj = @(x) At(h*x); % input = Ny*oy*Nx*ox vector; output = Ny*Nx matrix
fprintf('\nComputing covariance matrix...');
d = diag(h);
d = abs(d);
% d = ones(size(d)); % Disable weighting, simply do G^TPhi.
fprintf('\nPruning covariancemat ...');
nonzerocols = find(d >= prctile(d,100-klargestpercent));
d = d(nonzerocols);
d = max(diagthresholdepsilon, d);  % This ensures that inverting the values will not explode in computation
d12 = 1./sqrt(d);
Mask = sparse(1:length(nonzerocols), nonzerocols, ones(length(nonzerocols), 1), length(nonzerocols), (oy*Ny*ox*Nx));

% Final reduction operators Phi_grid and (Phi_grid)^T
B = @(x) d12.*(Mask*(h*A(x))); % subsampled, masked G^TPhi; input = [Ny, Nx] image; output = oy*Ny*ox*Nx vector
Bt = @(x) At((h*(Mask'*(d12.*x))));

evl = op_norm(B, Bt, [Ny, Nx], 1e-4, 200, verbosity);

% Embed the y using the same reduction
yTmat = d12.*(Mask*Gw'*y{1}{1});

epsilon = step_epsilon; % set epsilon value BEFORE running this script
epsilons = mat2cell(1.01*epsilon, length(epsilon));
epsilon = mat2cell(epsilon, length(epsilon));
epsilonT{1} = epsilon;
epsilonTs{1} = epsilons;

T = mat2cell([1], 1);
W = mat2cell(true(size(yTmat)), length(yTmat));
yT{1} = mat2cell(yTmat, length(yTmat));

%% PDFB parameter structure sent to the algorithm
param_pdfb.im = im; % original image, used to compute the SNR
param_pdfb.verbose = verbosity; % print log or not
param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb.nu2 = evl; % bound on the norm of the operator A*G
param_pdfb.gamma = 1e-5; % convergence parameter L1 (soft th parameter)
param_pdfb.tau = 0.49; % forward descent step size
param_pdfb.rel_obj = 1e-5; % stopping criterion
param_pdfb.max_iter = 10000; % max number of iterations
param_pdfb.lambda0 = 1; % relaxation step for primal update
param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb.sol_steps = [inf]; % saves images at the given iterations

param_pdfb.use_proj_elipse_fb = 1;
param_pdfb.elipse_proj_max_iter = 10;
param_pdfb.elipse_proj_min_iter = 1;
param_pdfb.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb.use_reweight_steps = 0;
param_pdfb.use_reweight_eps = 0;
param_pdfb.reweight_steps = [600:400:10000 inf];
param_pdfb.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb.reweight_min_steps_rel_obj = 50;
param_pdfb.reweight_alpha = 1; % Alpha always 1
param_pdfb.reweight_alpha_ff = 0.75; % 0.25 Too agressively reduces the weights, try 0.7, 0.8
param_pdfb.reweight_abs_of_max = inf;
param_pdfb.total_reweights = 20;

param_pdfb.use_adapt_bound_eps = 0;
param_pdfb.adapt_bound_steps = 100;
param_pdfb.adapt_bound_rel_obj = 1e-5;
param_pdfb.hard_thres = 0;
param_pdfb.adapt_bound_tol =1e-3;
param_pdfb.adapt_bound_start = 1000;

param_pdfb.savepath = save_path;

%% compute the solution
fprintf('Starting algorithms:\n\n');
tstart = tic;

script_run_all_tests_serial;

tend = toc(tstart);
fprintf('All algorithms runtime: %ds\n\n', ceil(tend));

%% save result and residual image
datasize = size(yT{1}{1}, 1);
imagesize = param_real_data.image_size_Ny * param_real_data.image_size_Nx;
fprintf('Data size:  %d\n', datasize);
fprintf('Image size: %d x %d\n', param_real_data.image_size_Ny, param_real_data.image_size_Nx);
fprintf('Data/Image = % 3d %%\n', round(100*(datasize/imagesize)));
filename = sprintf('results/rgrid/i%d.dl%d.k%03d.n%05d.eps%d.gamma%d.hyperspec.q6-5-51.rgrid.', param_real_data.image_size_Ny, (10*param_real_data.pixel_size), klargestpercent, param_pdfb.max_iter, step_epsilon, abs(log10(param_pdfb.gamma)));
predvisA = Gw*A(result_st.sol{1});
resvisA = y{1}{1} - predvisA;
resimgA = real(At(Gw'*resvisA));
delta = zeros(param_real_data.image_size_Ny, param_real_data.image_size_Nx); % create a delta function to compute the PSF
delta(param_real_data.image_size_Ny/2, param_real_data.image_size_Nx/2) = 1;
peakpsfA = max(max((real(At(h * A(delta)))))); % Peak of PSF = Peak(Phi^T(Phi(delta)))
resimgA = resimgA./peakpsfA;
fitswrite(fliplr(resimgA), sprintf('%sresA.fits', filename));
fitswrite(fliplr(result_st.sol{1}), sprintf('%srec.fits', filename));
predvisB =  B(result_st.sol{1});
resvisB = yT{1}{1} - predvisB;
resimgB = real(Bt(resvisB));
peakpsfB = max(max((real(Bt(B(delta)))))); % Peak of PSF = Peak(Phi'^T(Phi'(delta)))
resimgB = resimgB./peakpsfB;
fitswrite(fliplr(resimgB), sprintf('%sresB.fits', filename));
