
close all;

% file_name = 'M31s.fits';
% file_name = 'M31.fits';
% file_name = 'W28.fits';
file_name = 'galaxyCluster.fits';
% file_name = 'CYGCBEST.fits';
% file_name = 'cluster.fits';

NNy = 2048;
NNx = 2048;
Nfactor = 400;

Nps = 100;

range_sigmax = [2, 100];
range_sigmay = [2, 100];
range_sigma_cov = [0.1, 1.99];
range_amp = [1e-2 1];

prob_haze = 0.1;
range_haze_sigmax = [50, 120];
range_haze_sigmay = [50, 120];
range_haze_amp = [1e-16 1e-4];
range_haze_comp = [1 6];


[im, N, Ny, Nx] = util_read_image(file_name);

new_im = padarray(im, round([NNy - Ny, NNx - Nx]/2), 'both');

ps = zeros(NNy+Nfactor, NNx+Nfactor);


ps_posy = randi(NNy, Nps);
ps_posx = randi(NNx, Nps);
ps_sigmax = range_sigmax(1) + (range_sigmax(2) - range_sigmax(1)) * rand(Nps, 1);
ps_sigmay = range_sigmay(1) + (range_sigmay(2) - range_sigmay(1)) * rand(Nps, 1);
ps_sigma_cov = range_sigma_cov(1) + (range_sigma_cov(2) - range_sigma_cov(1)) * rand(Nps, 1);
ps_amp = range_amp(1) + (range_amp(2) - range_amp(1)) * rand(Nps, 1);


for i = 1:Nps
    mu = [ps_posy(i) ps_posx(i)];
    sigma = [ps_sigmay(i) ps_sigma_cov(i); ps_sigma_cov(i) ps_sigmax(i)];
    y = mu(1) - Nfactor:mu(1) + Nfactor;
    x = mu(2) - Nfactor:mu(2) + Nfactor;
    [Y,X] = meshgrid(y,x);

    F = mvnpdf([Y(:) X(:)], mu, sigma);
    F = reshape(F, length(y),length(x));
    F = ps_amp(i) * F/max(F(:));
    
    new_im(y(y>0 & y < NNy), x(x>0 & x < NNx)) = new_im(y(y>0 & y < NNy), x(x>0 & x < NNx)) + F(y>0 & y < NNy, x>0 & x < NNx);
end

figure; imagesc(log(new_im+1e-4));




