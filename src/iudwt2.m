function x = iudwt2(coef,W)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

N1 = length(W.dec);

t = 0;

for k = 1:N1
    [Nz Nx] = size(W.dec{k});
    t1 = Nz*Nx;
    W.dec{k}=reshape(coef(t+1:t+t1),Nz,Nx);
    t = t + t1;
end

x = indwt2(W);
end

