function [y0, y, y0f, yf, aY, input_snr_v, im_C] = util_gen_input_data(im, T, W, A, input_snr, use_per_block, snr_delta, uvidx, param)
% generates the input data


R = length(T);

if ~param.use_image_variation
    y0freq = A(im);
    im_C = class_ImageRI;
    im_C = im_C.setIm(im);
    im_C = im_C.setParams(param);
    
    y0 = cell(R, 1);
    for q = 1:R
        y0{q} = T{q} * y0freq(W{q});
    end
end

if param.use_image_variation
    if R ~= 1
        error('Multiple blocks are not supoerted when using a variable image \n\n');
    end
    
    for q = 1:R
        y0 = cell(R, 1);
        im_C = class_ImageRI;
        im_C = im_C.setIm(im);
        im_C = im_C.setParams(param);

        for k = 1:size(T{q}, 1)

            im_ = im_C.getIm(k);
            y0freq = A(im_);
            y0{q}(k) = T{q}(k, :) * y0freq(W{q});
            if mod(k, 5000) == 1
                fprintf('\nComputed %d out of %d visibility points with variable sky', k, size(T{q}, 1));
            end
        end
        fprintf('\nComputed %d out of %d visibility points with variable sky \n', size(T{q}, 1), size(T{q}, 1));
        y0{q} = y0{q}(:);
    end
end

Nm = numel(cell2mat(y0));
normy0 = norm(cell2mat(y0));

y = cell(R, 1);
aY = cell(R, 1);
input_snr_v = zeros(R, 1);

for q = 1:R
    Nmq = length(y0{q});
    
    if use_per_block
        input_snr_v(q) = input_snr + 2 * (rand(1, 1) - 0.5) * snr_delta;
    else
        input_snr_v(q) = input_snr;
    end
    
    % add Gaussian i.i.d. noise
    sigma_noise = 10^(-input_snr_v(q)/20) * normy0/sqrt(Nm);
    
    noise = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise/sqrt(2);
    y{q} = y0{q};
    y{q} = y{q} + noise;
    aY{q} = abs(y{q})./sigma_noise;
end
yf = cell2mat(y);
y0f = cell2mat(y0);

yf(cell2mat(uvidx)) = yf;
y0f(cell2mat(uvidx)) = y0f;

end





