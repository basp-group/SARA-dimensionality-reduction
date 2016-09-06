function im_out = grid2imgsize_fwd(im, st, ht, M)

% The forward operator implemented here is:
% Phi = S'Z'F'G'GFZS
% S = Scaling done to nullify the effect of the convolution by the
% interpolation kernel in G
% Z = Zero padding done before applying the FFT
% F = FFT
% G = Interpolation operator, maps grid Fourier points to continuous
% Fourier points (aka UV points). Each row of G has an interpolation kernel

% Input argument is an image of size [Ny, Nx]



% The structure st is different for the NUFFT toolbox used to have the
% Kaiser-Bessel kernel in the interpolation matrix
Ny1=st.Nd(1);
Nx1=st.Nd(2);
Ny2=st.Kd(1);
Nx2=st.Kd(2);

% Zero padding
imc = zeros(Ny2,Nx2);

% Find image first point: Note: No centered image while using the NUFFT package
% Zero padding done around the image, and image is not moved to the centre,
% because that's how the nufft_init function does it while simulating v0
% xo = floor(Nx2/2) - floor(Nx1/2);
% yo = floor(Ny2/2) - floor(Ny1/2);
xo = 0; yo = 0;
imc(yo+1:yo+Ny1,xo+1:xo+Nx1) = st.sn .* im; % Scaling

% FFT
spec = fft2(imc);%/sqrt(Nx2*Ny2); Multiplying factor removed, because in the gaussfwd_kbkern function 
% the normalisation of the FFT has been removed, 
% since the nufft_init function that we use to simulate v0 does not use the normalised FFT

% NO Quadrant swap after computing the FFT
% spec1 = fftshift(spec);

% Interpolation and gridding performed together with the holographic matrix
% (h = G'*G)
% M is a mask that exists if we trim/prune h to be smaller than the original size
% If h is not pruned then M = 1;
protospec = ht'*(M*spec(:)); % protospec lies on the gridded Fourier plane
% Inverse FFT
protoim = ifft2(reshape(protospec, Ny2, Nx2)); % protoim is on the (pseudo)image plane
% Cropping
% Find image first point
xo = 0; yo = 0;
im_out = conj(st.sn) .* protoim(yo+1:yo+Ny1,xo+1:xo+Nx1); % Cropping and scaling done in one step
im_out = im_out(:);
