function mask = genmask(pdf, seed)
% GENMASK - Generate a mask with variable density sampling
%
% Inputs:
% pdf : Sampling profile.
% seed : Seed for the random number generator. Optional.

if nargin==2
    rand('seed', seed);
end

mask = rand(size(pdf))<pdf;

