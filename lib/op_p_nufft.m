function [A, At, G, W, Gw, As, Ats, S] = op_p_nufft(p, N, Nn, No, Ns, ww, param)

% Create the nonuniform gridding matrix and fft operators to be used for
% parallel processing
%
% in:
% p{:}[2] - nonuniformly distributed frequency location points for each
%           cell member which will be treated in parallel
% N[2]    - size of the reconstruction image
% Nn[2]   - size of the kernels (number of neighbors considered on each direction)
% No[2]   - oversampled fft from which to recover the non uniform fft via
%           kernel convolution
% Ns[2]   - fft shift
%
% out:
% A[@]          - function handle for direct operator
% At[@]         - function handle for adjoint operator
% G{:}[:][:]    - convolution kernel matrix (small) associated with each
%               patch in the fourier plane
% W{:}          - mask of the values that contribute to the convolution
% Gw[:][:]      - global convolution kernel matrix

if ~exist('param', 'var')
    param = struct();
end
if ~isfield(param, 'use_nufft_blocks'), param.use_nufft_blocks = 1; end
if ~isfield(param, 'gen_only_fft_op'), param.gen_only_fft_op = 0; end
if ~isfield(param, 'gen_fft_op_without_scale'), param.gen_fft_op_without_scale = 0; end
if ~exist('ww', 'var')
    ww = cell(length(p), 1);
    for q=1:length(p)
        ww{q} = ones(length(p{q}(:, 1)), 1);
    end
end
if ~isfield(param, 'use_fft_mask'), param.use_fft_mask = 1; end

R = size(p, 1);

if param.gen_fft_op_without_scale
    [As, Ats, ~, S] = op_nufft([0, 0], N, Nn, No, Ns, param.gen_fft_op_without_scale);
else
    As = [];
    Ats = [];
    S = [];
end

if param.gen_only_fft_op
    [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns, 0);
    G = [];
    W = [];
    Gw = [];
else
    if ~param.use_nufft_blocks
        %% compute the overall gridding matrix and its associated kernels
        [A, At, Gw, ~] = op_nufft(cell2mat(p), N, Nn, No, Ns, 0);

        %% compute small gridding matrices associated with each parallel block
        G = cell(R, 1);
        if param.use_fft_mask
            W = cell(R, 1);
        else
            W = cell(R, 1);
            for q = 1:R
                W{q} = ':';
            end
        end

        % block start position
        fprintf('\nComputing block matrices ...\n');
        b_st = 1;
        for q = 1:R
            tstart = tic;
            % current block length
            % the matrix Gw is structured identical to the structure of p thus we 
            % grab it block by block
            b_l = length(p{q});

            % get a block out of the large G and trim it
            Gw(b_st:b_st+b_l-1, :) = spdiags(ww{q}, 0, b_l, b_l) * Gw(b_st:b_st+b_l-1, :);
            Gb = Gw(b_st:b_st+b_l-1, :);

            %% now trim the zero rows and store a mask in W
            
            if param.use_fft_mask
                % preallocate W for speed
                W{q} = false(No(1)*No(2), 1);
            end

            if param.use_fft_mask
                % use the absolute values to speed up the search
                Gb_a = abs(Gb);
            
                % check if eack line is entirely zero
                W{q} = Gb_a' * ones(b_l, 1) ~= 0;
            end
            
            if param.use_fft_mask
                % store only what we need from G
                G{q} = Gb(:, W{q});
            else
                G{q} = Gb;
            end

            % iterate among the blocks
            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end
    else

        %% compute small gridding matrices associated with each parallel block
        
        Gw = spalloc(length(cell2mat(p)), No(1)*No(2), 16 * length(cell2mat(p)));
        G = cell(R, 1);
        W = cell(R, 1);

        b_st = 1;
        % block start position
        fprintf('\nComputing block matrices ...\n');
        for q = 1:R

            tstart = tic;
            b_l = length(p{q});


            %% compute the small gridding matrix and its associated kernels
            [~, ~, Gb, ~] = op_nufft([p{q, 1} p{q, 2}], N, Nn, No, Ns, 0);

            %% now trim the zero rows and store a mask in W

            if param.use_fft_mask
                % preallocate W for speed
                W{q} = false(No(1)*No(2), 1);
            end

            Gb = spdiags(ww{q}, 0, b_l, b_l) * Gb;
            
            if param.use_fft_mask
                % use the absolute values to speed up the search
                Gb_a = abs(Gb);

                % check if eack line is entirely zero
                W{q} = Gb_a' * ones(size(Gb, 1), 1) ~= 0;
            end

            if param.use_fft_mask
                % store only what we need from G
                G{q} = Gb(:, W{q});
            else
                G{q} = Gb;
            end

            %% fill the whole Gw
            Gw(b_st:b_st+b_l-1, :) = Gb;

            b_st = b_st+b_l;
            tend = toc(tstart);
            fprintf('Block matrix %d: %ds \n', q, ceil(tend));
        end

        [A, At, ~, ~] = op_nufft([0, 0], N, Nn, No, Ns, 0);
    end
end


end
