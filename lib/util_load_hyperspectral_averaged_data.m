function [y, uw, vw, sigma_noise, nW] = util_load_hyperspectral_averaged_data(dirname, param)

if ~isfield(param, 'frequencyaveraging'), param.frequencyaveraging = 0; end

startlimit = 0; endlimit = 63; % dirname contains 64 channels, numbered 0 to 63
% Shouldn't run beyond the limit of channel indices

y = []; uw = []; vw = []; nW = [];
u = []; v = [];

% The hyperspectral data is read in the following way:
% 1. Identify the central frequency (middle of the band)
% 2. See how many channels we want to squish together to make the visibilities dataset
% 3. Check the desired gap between channels, if any
% 4. Average the visibilities at a UV point
% 5. Scale the uv points to fit the [-pi, pi] range, taking care of maximum baselines
numfrequencies = 10; % actually we have 64 channels in data/vis/WEIGHTED-CYGA-C-6680-64CH
centralfrequency = 31;
frequencystep = 5;
startfrequency = centralfrequency - (frequencystep * (floor(numfrequencies/2)));
endfrequency = centralfrequency + (frequencystep * (ceil(numfrequencies/2) - 1));
if ((startfrequency < startlimit) || (endfrequency > endlimit))
    error('Incompatible frequency steps and channel numbers in %s', dirname);
end
for freq = startfrequency:frequencystep:endfrequency % skip channels - taking all would be too much data together at the moment
    % freq = 0;
    currfreq = sprintf('freq%d', freq);
    currfilename = sprintf('%s/weighted_I_%d.mat', dirname, freq);
    currdata = load(currfilename);

    uvw.(currfreq) = currdata.uvw;
    y_I.(currfreq) = (currdata.y_I)';
    sigma.(currfreq) = (currdata.sigmas)';
    flag.(currfreq) = (currdata.flag)';
    weights.(currfreq) = (currdata.weights)';

    % Remove flagged entries and entries with sigma = 0
    usable_indices = (flag.(currfreq) == 0) & (sigma.(currfreq) > 1e-4);
    y_usable = double(y_I.(currfreq)(usable_indices));
    uvw_usable = uvw.(currfreq)(usable_indices, :);
    weights_usable = weights.(currfreq)(usable_indices);
    weights_usable = sqrt(weights_usable);

    sigma_noise = 1;
    if(~param.frequencyaveraging)
        y = [y ; (y_usable .* weights_usable)];
        u = [u ; uvw_usable(:,1)];
        v = [v ; uvw_usable(:,2)];
        nW = [nW ; weights_usable];
    else
        if freq == startfrequency
            y = [y ; (y_usable .* weights_usable)];
        else
            y = y + (y_usable .* weights_usable);
        end
        if freq == centralfrequency
            u = [u ; uvw_usable(:,1)];
            v = [v ; uvw_usable(:,2)];
            nW = [nW ; weights_usable];
        end
    end
end
if(param.frequencyaveraging)
    y = y/numfrequencies;
end
bmaxProj = max(sqrt(u.^2+v.^2));
dl = param.pixel_size;
vw = v * pi / (bmaxProj * dl);
uw = u * pi / (bmaxProj * dl);


