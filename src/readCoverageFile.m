function [u, v, w] = readCoverageFile(coveragefile, numUVpoints)

%readCoverageFile

%addpath ../data % Because coverage files are usually stored in the data folder
u = []; v = []; w = [];


% Bool: If the UV points were already generated and stored in a MAT file.
% To be careful that the number of UV points in the MAT file matches what
% we want
saveduvpoints = 1;

% This was my implementation, simply reading from the file line-by-line
%     % Load uv coverage from the Meqtree simulated file obtained from Arwa
%     meerkat = readtable(coveragefile,'Delimiter',' ','ReadVariableNames',true,'Format','%f%f%f');
%     u = meerkat{:,1};
%     maxabsu = max(abs(max(u)), abs(min(u)));
%     u = u*pi/maxabsu;
%     v = meerkat{:,2};
%     maxabsv = max(abs(max(v)), abs(min(v)));
%     v = v*pi/maxabsv;
    
if(saveduvpoints)
    load(coveragefile);
    u = uvw(:,1); % When a mat file generated from Pyxis is loaded,
    v = uvw(:,2); % the variable is uvw, with three columns
    w = uvw(:,3);
else
% This is Arwa's implementation, with more info about the simulation
    [u, v, w, R] = util_uv_read(coveragefile, numUVpoints/2); % Divide by 2 because in the function the UV points are generated as mirror images
end

%assert(length(u) >= numUVpoints, 'Number of UV points asked: %d\nNumber of UV points loaded from file: %d', numUVpoints, length(u));

end

