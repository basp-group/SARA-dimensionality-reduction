function [x] = vec2complex(v)
    l = length(v);
    x = v(1:l/2) + 1i*v(l/2+1:end);
end