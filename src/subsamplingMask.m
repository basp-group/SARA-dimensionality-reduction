function [fatmatrix] = subsamplingMask(reducedSize, origSize)
% subsamplingMask: Create a fat matrix which randomly subsamples
% from a long vector to give a short vector
% We could very well use datasample() for this, but we need an explicit
% instance of a matrix, which we will transpose and use in other places
    nonzerocols = sort(randperm(origSize, reducedSize));
    fatmatrix = sparse(1:reducedSize, nonzerocols, ones(reducedSize, 1), reducedSize, origSize);
end

