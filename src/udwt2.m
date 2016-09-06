function coef = udwt2(x,W)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


W1 = ndwt2(x,W.level,W.filters,'mode',W.mode);

N1=length(W1.dec);
coef = [];

for i = 1:N1
    coef = [coef; W1.dec{i}(:)];
end

end

