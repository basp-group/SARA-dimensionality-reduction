function [A, At] = op_nu_so_fft2(N, No, scale, use_fft_on_gpu)
% Oversampled ftt2 and scaled operator computed by a modified nufft function
%
% in:
% N[2]     - size of the reconstruction image
% No[2]    - oversampled fft from which to recover the non uniform fft via
% scale[:] - scale parameters precomputed by nufft_init
%
% out:
% A[@]     - function handle for direct operator
% At[@]    - function handle for adjoint operator

if ~exist('use_fft_on_gpu', 'var')
    use_fft_on_gpu = 0;
end

if use_fft_on_gpu == 0
    A = @(x) so_fft2(x, No, scale);
    At = @(x) so_fft2_adj(x, N, No, scale);
else
    A = @(x) so_fft2_gpu(x, No, scale);
    At = @(x) so_fft2_adj_gpu(x, N, No, scale);
end

end


