close all
%%
param.use_shift = 0;
param.pixel_size = 1;
param.use_undersamplig = 0;

Ny = 2^7;
Nx = 2^7;

% data = 'data/3C129-vla-full.mat';
% data = 'data/vla-3C129-luke-test.mat';
data = 'data/vla-3C129-luke-pruned-500-noise.mat';

[y, uw, vw, ~, ~] = util_load_real_vis_data(data, param);

%%
[aY, aUV, aPos] = util_gen_uv_dist_per_pixel(uw, vw, Ny, Nx, y);

%%
IMr = zeros(Ny, Nx);
for k=1:Ny
    for j=1:Nx
%         value = -min(real(aY{k, j})) + max(real(aY{k, j}));
%         value = std(abs(real(aY{k, j})))/abs(mean(abs(real(aY{k, j}))));
%         value = std((real(aY{k, j})))/abs(mean((real(aY{k, j}))));
        value = mean(real(aY{k, j}));
        if isempty(value)
            value = 0;
        end
        IMr(k, j) = value;
    end
end
figure(1); imagesc(IMr); colorbar;


IMi = zeros(Ny, Nx);
for k=1:Ny
    for j=1:Nx
%         value = -min(imag(aY{k, j})) + max(imag(aY{k, j}));
%         value = std(abs(imag(aY{k, j})))/abs(mean(abs(imag(aY{k, j}))));
%         value = std((imag(aY{k, j})))/abs(mean((imag(aY{k, j}))));
        value = mean(imag(aY{k, j}));
        if isempty(value)
            value = 0;
        end
        IMi(k, j) = value;
    end
end
figure(2); imagesc(IMi); colorbar;



%% 
% remove points
th = 0.6;

figure(3);
mr = max(IMr(:));


allPos = zeros(length(y), 1);

st = 0;

for k=1:Ny
    for j=1:Nx
        if IMr(k, j) > th * mr
            scatter(aUV{k, j}(:, 1), aUV{k, j}(:, 2), '.'); hold on;
            allPos(st+1:st+length(aPos{k, j})) = aPos{k, j};
            st = st + length(aPos{k, j});
        end
    end
end

th = 0.6;

figure(4);
mi = max(IMi(:));

for k=1:Ny
    for j=1:Nx
        if IMi(k, j) > th * mi
            scatter(aUV{k, j}(:, 1), aUV{k, j}(:, 2), '.'); hold on;
            allPos(st+1:st+length(aPos{k, j})) = aPos{k, j};
            st = st + length(aPos{k, j});
        end
    end
end

allPos(st:end) = [];
allPos = unique(allPos);

%%
if false
    %%
    load(data);
    pos = allPos;
    uvw(pos, :) = [];
    weights(pos) = [];
    y_I(pos) = [];
    y_V(pos) = [];
    save('data/vla-3C129-luke-pruned-500-noise-rem.mat', 'uvw', 'weights', 'y_I', 'y_V');
end

