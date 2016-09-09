% Test file to to run image reconstruction with embedded data
% Dimensionality reduction performed with Rsing and Rgrid as defined in
% "A Fourier dimensionality reduction model for big data in interferometric imaging" by S. Vijay Kartik et al.
% Author: S. Vijay Kartik


% 256x256 image of M31 (N = 256*256)
imgfile = 'M31_256.fits';
% Coverage loaded from file
coveragefile = 'ska254.i256.p10.uvw.mat';
% Initial data size of 10N
visibSize = 256*256*10;
% Input SNR (in dB) to calculate additive i.i.d. Gaussian noise
input_snr = 30;
% For further below-image size discarding for Rsing and Rgrid, 
% retain only largest K percent of values
klargestpercent = 5;
% Run number for current simulation
run = 1;

[SNR_allvisibs, DR_allvisibs, M_allvisibs, ~, tend_allvisibs] = allvisibs_admm(visibSize, input_snr, imgfile, coveragefile, run);

[SNR_rsing, DR_rsing, M_rsing, ~, tend_rsing] = rsing_admm(visibSize, input_snr, imgfile, coveragefile, run);
[SNR_rsing_discarding, DR_rsing_discarding, M_rsing_discarding, ~, tend_rsing_discarding] = rsing_discardsingvalues_admm(visibSize, input_snr, imgfile, coveragefile, klargestpercent, run);

[SNR_rgrid, DR_rgrid, M_rgrid, ~, tend_rgrid] = rgrid_admm(visibSize, input_snr, imgfile, coveragefile, run);
[SNR_rgrid_discarding, DR_rgrid_discarding, M_rgrid_discarding, ~, tend_rgrid_discarding] = rgrid_discardsingvalues_admm(visibSize, input_snr, imgfile, coveragefile, klargestpercent, run);


