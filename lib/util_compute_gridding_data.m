function [yT, T, Tw] = util_compute_gridding_data(y, G, Gw)
% Create the nonuniform holographic matrix and fft operators to be used for
% parallel processing. G'G case 
%
% in:
% G{:}[:][:]    - convolution kernel matrix
% Gw[:][:]      - global convolution kernel matrix
%
% out:
% T{:}[:][:]    - gridded convolution kernel matrix
% Tw[:][:]      - gridded global convolution kernel matrix

%% compute small gridding matrices associated with each parallel block G'G
R = length(G);
num_tests = length(y);
yT = cell(num_tests, 1);
T = cell(R, 1);

fprintf('\nComputing block square matrices G''G ... \n');

tstart = tic;
Tw = Gw' * Gw;
tend = toc(tstart);
fprintf('Whole G''G: time %ds \n', ceil(tend));


for q = 1:R
    tstart = tic;
    T{q} = G{q}' * G{q};
    tend = toc(tstart);
    fprintf('Square block matrix %d: %ds \n', q, ceil(tend));
    
    for k = 1:num_tests
        yT{k}{q} = G{q}' * y{k}{q};
    end
end
for k = 1:num_tests
    yT{k} = yT{k}';
end

end

