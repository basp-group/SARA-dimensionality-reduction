function [u,v] = simulateCoverage(numPixels, numUVpoints)
%simulateCoverage

% Boolean to know which coverage to simulate
isholeyuniformcoverage = 1;
isholeygaussiancoverage = 0;
isholeygeometriccoverage = 0;
% Simulate coverage using a Gaussian sampling with random holes.
% All parameters of hole sizes/numbers and others are set here as well

if(isholeygaussiancoverage)
    util_gen_sampling_pattern_config; % Set all parameters
    sparam.N = numPixels;
    sparam.p = ceil(numUVpoints/numPixels);
    sampling_pattern = 'gaussian+large-holes';
    [u, v, uw, vw, Nm] = util_gen_sampling_pattern(sampling_pattern, sparam);
    u = u{1}; v = v{1};
elseif(isholeyuniformcoverage)
    util_gen_sampling_pattern_config; % Set all parameters
    sparam.N = numPixels;
    sparam.p = ceil(numUVpoints/numPixels);
    sampling_pattern = 'uniform+large-holes';
    [u, v, uw, vw, Nm] = util_gen_sampling_pattern(sampling_pattern, sparam);
    u = u{1}; v = v{1};
elseif(isholeygeometriccoverage)
    util_gen_sampling_pattern_config; % Set all parameters
    sparam.N = numPixels;
    sparam.p = ceil(numUVpoints/numPixels);
    sampling_pattern = 'geometric+large-holes';
    [u, v, uw, vw, Nm] = util_gen_sampling_pattern(sampling_pattern, sparam);
    u = u{1}; v = v{1};
else
% Simulate coverage using a Gaussian sampling pattern on the uv plane
% within [-pi, pi]
% Sampling pattern
    sigma_m = pi/6;
    rho = 2-(erf(pi/(sigma_m*sqrt(2))))^2;
    num_meas1 = floor(numUVpoints*rho);
    u1 = sigma_m*randn(num_meas1,1);
    v1 = sigma_m*randn(num_meas1,1);
    %Discard points outside (-pi,pi)x(-pi,pi)
    sf1=find((u1<pi)&(u1>-pi));
    sf2=find((v1<pi)&(v1>-pi));
    sf=intersect(sf1,sf2);
    v=v1(sf);
    u=u1(sf);
    clear v1 u1
end

end
