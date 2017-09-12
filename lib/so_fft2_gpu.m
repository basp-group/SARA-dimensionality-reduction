function X = so_fft2_gpu(x, No, scale)
% Computes and oversampled and scaled FFT2 with the scale factors
% precomputed from nufft_init
%
% in:
% x[:][:]  - input image
% No[2]    - overscale size
% scale[:] - scale parameters precomputed by nufft_init
%
% out:
% X[:]     - FFT2 coefficients


% apply scaling factors
x = x .* scale; 

% oversampled FFT 
X = gather(fft2(gpuArray(x), No(1), No(2)));
X = X(:);

end