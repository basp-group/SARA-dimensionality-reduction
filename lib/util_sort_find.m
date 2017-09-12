function [lower_index,upper_index] = util_sort_find(x,LowerBound,UpperBound,exactLower,exactUpper)
% fast O(log2(N)) computation of the range of indices of x that satify the
% upper and lower bound values using the fact that the x vector is sorted
% from low to high values. Computation is done via a binary search.
%
% Input:
%
% x-            A vector of sorted values from low to high.       
%
% LowerBound-   Lower boundary on the values of x in the search
%
% UpperBound-   Upper boundary on the values of x in the search
%
% Output:
%
% lower_index-  The smallest index such that
%               LowerBound<=x(index)<=UpperBound
%
% upper_index-  The largest index such that
%               LowerBound<=x(index)<=UpperBound

x = [-inf; x(:); inf];

if LowerBound>x(end) || UpperBound<x(1) || UpperBound<LowerBound
    % no indices satify bounding conditions
    lower_index = [];
    upper_index = [];
    return;
end

lower_index_a=1;
lower_index_b=length(x); % x(lower_index_b) will always satisfy lowerbound
upper_index_a=1;         % x(upper_index_a) will always satisfy upperbound
upper_index_b=length(x);

%
% The following loop increases _a and decreases _b until they differ 
% by at most 1. Because one of these index variables always satisfies the 
% appropriate bound, this means the loop will terminate with either 
% lower_index_a or lower_index_b having the minimum possible index that 
% satifies the lower bound, and either upper_index_a or upper_index_b 
% having the largest possible index that satisfies the upper bound. 
%
while (lower_index_a+1<lower_index_b) || (upper_index_a+1<upper_index_b)

    lw=floor((lower_index_a+lower_index_b)/2); % split the upper index

    if x(lw) >= LowerBound
        lower_index_b=lw; % decrease lower_index_b (whose x value remains \geq to lower bound)   
    else
        lower_index_a=lw; % increase lower_index_a (whose x value remains less than lower bound)
        if (lw>upper_index_a) && (lw<upper_index_b)
            upper_index_a=lw;% increase upper_index_a (whose x value remains less than lower bound and thus upper bound)
        end
    end

    up=ceil((upper_index_a+upper_index_b)/2);% split the lower index
    if x(up) <= UpperBound
        upper_index_a=up; % increase upper_index_a (whose x value remains \leq to upper bound) 
    else
        upper_index_b=up; % decrease upper_index_b
        if (up<lower_index_b) && (up>lower_index_a)
            lower_index_b=up;%decrease lower_index_b (whose x value remains greater than upper bound and thus lower bound)
        end
    end
end

if x(lower_index_a)>=LowerBound
    lower_index = lower_index_b;
else
    lower_index = lower_index_a;
end
if x(upper_index_b)<=UpperBound
    upper_index = upper_index_a;
else
    upper_index = upper_index_b;
end

lower_index = lower_index + 1;
upper_index = upper_index - 1;

while ~exactLower && lower_index <= length(x) && x(lower_index) == LowerBound
    lower_index = lower_index + 1;
end

while ~exactUpper && upper_index >= 1 && x(upper_index) == UpperBound
    upper_index = upper_index - 1;
end

lower_index = lower_index - 1;
upper_index = upper_index - 1;