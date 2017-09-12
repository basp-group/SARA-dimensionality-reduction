function [x,xa] = util_rand_2dggd(c,s,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates complex generalized gaussian random variables with augmented
% covariance matrix Ta = [2*s 0; 0 2*s];
% and shape parameter c, where c = 1 corresponds to the Gaussian case.
% x = util_rand_2dggd(c,s,N) generates a vector 2xN of samples with a circular
% gaussian distribution with shape parameter c, and variance 2*s.
% [x,xa] = util_rand_2dggd(c,s,N) generates an additional augmented matrix 2xN,
% xa = [x;conj(x)].
% The samples are generated according to the results presented in:
% Mike Novey, T. Adali, and A. Roy "A complex generalized Gaussian Distribution---
% characterization, generation, and estimation" in
% IEEE Trans. Signal Proc., Vol 58. No 3, MARCH 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
c = c/2;
zn = (gamrnd(1/c,1,1,N).^(1/(2*c))).*(exp(sqrt(-1)*2*pi*rand(1,N)));
eta = gamma(2/c)/gamma(1/c);
Wn = zn./sqrt(eta);
W = [Wn;conj(Wn)];
Ta = sqrtm([2*s 0; 0 2*s]);
xa = Ta*W;
x = xa(1,:);
x = [real(x); imag(x)]';