function [SNR, DR, M, histogrampeakiness, tend] = allvisibs_admm(visibSize, input_snr, imgfile, coveragefile, run)
% inputs
% 1. visibSize = total number of visibilities in the data vector [usually some
% multiple of signal vector size p*(imgSize)]
% 2. input_snr = input signal to noise ratio, to define the sigma of the additive noise in the measurement
% 3. imgfile = pathname of the image file on disk
% 4. coveragefile (optional), if loading the coverage from disk. Usually the
% coverage is generated on the fly
% 5. run = run number to tag the simulation instance (log and result files are uniquely named with this run number)

addpath data
addpath data/images
addpath src
addpath src/irt
try
    setup;
catch ME
    error('NUFFT library not found in location src/irt');
%     rethrow(ME);
end

rng('shuffle');

%% Reading input image
% Moved the reading image routine to the very top for the following reasons
% 1. First check if image can be opened. If not, exit - no need to continue
% 2. Compute image size to give correct names to log and results files
fprintf('Reading image file... ');
im = readimg(imgfile); % This will exit with error if fitsread fails to read  the file
[Ny, Nx]=size(im);
N = Ny*Nx;
fprintf('Done\n');

%% Find appropriate directory, open log and results files
% Write log file and result file in the specific subdirectory depending on
% the test image
if regexpi(imgfile, '\w*m31\w*.fits')
    subdirname = 'm31';
elseif regexpi(imgfile, '\w*30dor\w*.fits')
    subdirname = '30dor';
elseif regexpi(imgfile, '\w*M87ABABCPCONV6\w*.fits')
    subdirname = 'm87';
elseif regexpi(imgfile, '\w*g41\w*.fits')
    subdirname = 'g41';
elseif regexpi(imgfile, '\w*CYG\w*.fits')
    subdirname = 'cyga';
elseif regexpi(imgfile, '\w*cluster\w*.fits')
    subdirname = 'galaxycluster';
else
    error('Unknown test image - no corresponding log/results folder');
end
% Filenames to save the solution, text output and relevant variables
solfile=sprintf('results/%s/soln.allvisibs.i%d.p%02d.r%02d.mat', subdirname, Ny, ceil(visibSize/(N)), run);
logfile=sprintf('logs/%s/soln.allvisibs.i%d.p%02d.r%02d.log', subdirname, Ny, ceil(visibSize/(N)), run);
diary(logfile);

% Flag to cheat in the computation of epsilon by getting the EXACT value thanks to
% the noise realisation used in the simulation (God-like powers)
cheatepsilon = 0;
% Flag to set if we want to plot the reconstructed, dirty and residual images
seeplots = 1;
% Flag to use simulated coverage or pre-existing 'real' telescope coverage
coveragefileexists = 1;
% Flag to use precomputed matrix G
Gfileexists = 0;

% Get uv coverage
tstart1=tic;
if(coveragefileexists)
    fprintf('Loading coverage file... ');
    [u1, v1, w1] = readCoverageFile(coveragefile, visibSize);
    zoomthreshold = 6e15;
    u = u1(abs(u1)<zoomthreshold & abs(v1)<zoomthreshold);
    v = v1(abs(u1)<zoomthreshold & abs(v1)<zoomthreshold);
    scaler = max(max(u), max(v));
    u = u*pi/scaler; v = v*pi/scaler;
else
    fprintf('Simulating coverage... ');
    [u,v] = simulateCoverage(N, visibSize);
end
tend1=toc(tstart1);
fprintf('Done\n');
fprintf('Time: %e\n', tend1);
%% Measurement operator initialization 

%Oversampling factors for nufft
ox = 2;
oy = 2;

%Number of neighbours for nufft
Kx = 8;
Ky = 8;

% Get degridding matrix
tstart1=tic;
if(Gfileexists)
    fprintf('Loading G matrix file... ');
    load(Gfile);
else
    %Initialize nufft parameters
    fprintf('Initializing the NUFFT operator... ');
    st = nufft_init([v u],[Ny Nx],[Ky Kx],[oy*Ny ox*Nx], [Ny/2 Nx/2]);
end
tend1=toc(tstart1);
fprintf('Done\n');
fprintf('Time: %e\n', tend1);

A = @(x) nufft(x, st);
At = @(x) nufft_adj(x, st);

%Maximum eigenvalue of operator A^TA
eval = pow_method(A, At, [Ny,Nx], 1e-4, 100, 1);

% Simulate noisy data
fprintf('Simulating noisy data... ');
y0 = A(im);
    
M = length(y0);
    
% Add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20)*norm(y0(:))/sqrt(M);
noise = (randn(size(y0)) + 1i*randn(size(y0)))*sigma_noise/sqrt(2);
y = y0 + noise;

fprintf('Done\n');

%% Parameter estimation

%Bound for the L2 norm
fprintf('Computing epsilon bound... ');
if(cheatepsilon)
    %%% compute epsilon from direct calculation with the noise realisation
    epsilon = norm(noise);
else
    epsilon = sqrt(M + 2*sqrt(M))*sigma_noise;
end

histogrampeakiness = 0; % dummy output argument, just to be consistent with the output of the D^(-1/2) cholesky approximation reconstruction
fprintf('Done\n');

%Dirty image
fprintf('Computing dirty image... ');
dirty = At(y);
dirty1 = 2*real(dirty)/eval;
fprintf('Done\n');

%% Sparsity operator definition

%Wavelets parameters
nlevel=4;
dwtmode('per');

% Sparsity operator for SARA

[C1,S1]=wavedec2(dirty1,nlevel,'db1'); ncoef1=length(C1);
[C2,S2]=wavedec2(dirty1,nlevel,'db2'); ncoef2=length(C2);
[C3,S3]=wavedec2(dirty1,nlevel,'db3'); ncoef3=length(C3);
[C4,S4]=wavedec2(dirty1,nlevel,'db4'); ncoef4=length(C4);
[C5,S5]=wavedec2(dirty1,nlevel,'db5'); ncoef5=length(C5);
[C6,S6]=wavedec2(dirty1,nlevel,'db6'); ncoef6=length(C6);
[C7,S7]=wavedec2(dirty1,nlevel,'db7'); ncoef7=length(C7);
[C8,S8]=wavedec2(dirty1,nlevel,'db8'); ncoef8=length(C8);

clear C1 C2 C3 C4 C5 C6 C7 C8


Psit = @(x) [wavedec2(x,nlevel,'db1')'; wavedec2(x,nlevel,'db2')';...
    wavedec2(x,nlevel,'db3')';wavedec2(x,nlevel,'db4')';...
    wavedec2(x,nlevel,'db5')'; wavedec2(x,nlevel,'db6')';...
    wavedec2(x,nlevel,'db7')'; wavedec2(x,nlevel,'db8')'; x(:)]/sqrt(9); 
Psi = @(x) (waverec2(x(1:ncoef1),S1,'db1')+...
    waverec2(x(ncoef1+1:ncoef1+ncoef2),S2,'db2')+...
    waverec2(x(2*ncoef1+1:2*ncoef1+ncoef2),S3,'db3')+...
    waverec2(x(3*ncoef1+1:3*ncoef1+ncoef2),S4,'db4')+...
    waverec2(x(4*ncoef1+1:4*ncoef1+ncoef2),S5,'db5')+...
    waverec2(x(5*ncoef1+1:5*ncoef1+ncoef2),S6,'db6')+...
    waverec2(x(6*ncoef1+1:6*ncoef1+ncoef2),S7,'db7')+...
    waverec2(x(7*ncoef1+1:7*ncoef1+ncoef2),S8,'db8')+...
    reshape(x(8*ncoef1+1:8*ncoef1+N), [Ny Nx]))/sqrt(9);


fprintf('Running solver...\n');
% Parameters for BPDN
param1.verbose = 1; % Print log or not
param1.gamma = 1e-6; % Converge parameter
param1.rel_obj = 1e-4; % Stopping criterion for the L1 problem
param1.max_iter = 100000; % Max. number of iterations for the L1 problem
param1.nu = eval; % Bound on the norm of the operator A
param1.tight_L1 = 0; % Indicate if Psit is a tight frame (1) or not (0)
param1.max_iter_L1 = 100;
param1.rel_obj_L1 = 1e-2;
param1.pos_L1 = 1;
param1.nu_L1 = 1;
param1.verbose_L1 = 0; % Print log or not

     
%Solve BPDN
tstart2 = tic;
[sol, z] = admm_bpcon(y, epsilon, A, At, Psi, Psit, param1);
tend = toc(tstart2);

err = im - sol;
SNR = 20*log10(norm(im(:))/norm(err(:)));

residual = At(y - A(sol));
DR = eval*max(sol(:))/(norm(residual(:))/sqrt(N));

fprintf('SNR: %f\n', SNR);
fprintf('DR: %f\n', DR);
fprintf('Visibilities dimension: %d\n', M);
fprintf('Reconstruction time: %f seconds\n', tend);

soln = sol/max(sol(:));
dirtyn = dirty1 - min(dirty1(:));
dirtyn = dirtyn/max(dirtyn(:));

fprintf('Saving reconstructed, residual and dirty images... ');
save(solfile, 'soln', 'residual', 'dirtyn', 'eval');
fprintf('Done\n');

if(seeplots && ismac)
    figure, imagesc(log10(soln + 1e-4)), colorbar, axis image, title('Reconstructed');
    figure, imagesc(log10(dirtyn + 1e-3)), colorbar, axis image, title('Dirty');
    figure, imagesc(real(residual)/eval), colorbar, axis image, title('Residual');
end
