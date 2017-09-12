

    if gen_figures
        fprintf('   Generating figures and computing solution sparsity ... \n');
        snr = 0;
        sp = 0;
        soln = zeros(Ny, Nx);
        residualn = zeros(Ny, Nx);
        dirtyn = zeros(Ny, Nx);
        for i = 1:num_tests

            asol = result_st.sol{i};
            ay = y{i};
            aL1_v = result_st.L1_v{i};
            aL2_v = result_st.L2_v{i};
            aL1_vp = result_st.L1_vp{i};
            aL2_vp = result_st.L2_vp{i};
%             ay0 = y0{i};
            aepsilon = epsilon{i};
            aepsilonT = epsilonT{i};
            
            if ~use_real_visibilities
                error = im - asol;
                snr_ = 20 * log10(norm(im(:))/norm(error(:)));
                snr = snr + snr_;
            end

            R = length(ay);
            ys = A(asol);



            if exist('Gw', 'var')
                residual = At(Gw' * (cell2mat(ay) - Gw * ys));
            end
            if ~exist('Gw', 'var') && exist('G', 'var')
                residual = zeros(size(W{1}));
                for q = 1:R
                    residual(W{q}) = residual(W{q}) + G{q}' * (ay{q} - G{q} * ys(W{q}));
                end
                residual = At(residual);
            end
            residualn_ = real(residual)/evl;
            residualn = residualn + residualn_;

            soln_ = asol;%/max(asol(:));
            soln = soln + soln_;
            if ~use_real_visibilities
                imn = im/max(im(:));
            end

            wcoef = [];
            for k = 1:length(Psit)
                wcoef = [wcoef; Psit{k}(asol)];
            end

            sp_ = sum(abs(wcoef) > 1e-3 * max(abs(wcoef)));
            sp = sp + sp_;

            No = length(ys);

            % dirty image
            if exist('Gw', 'var')
                dirty = At(Gw' * cell2mat(ay));
            end
            if ~exist('Gw', 'var') && exist('G', 'var')
                dirty = zeros(size(W{R}));
                for q = 1:R
                    dirty(W{q}) = dirty(W{q}) + G{q}' * ay{q};
                end
                dirty = At(dirty);
            end


            dirty = 2 * real(dirty);

            % normalized dirty image
            dirtyn_ = dirty - min(dirty(:));
            dirtyn_ = dirtyn_/max(dirtyn_(:));
            dirtyn = dirtyn + dirtyn_;



            if ~gen_only_average_figures

                fprintf('\nTest: %i \n', i);
                if ~use_real_visibilities
                    fprintf('\n\n    SNR: %f \n', snr_);
                end
                fprintf('    Solution sparsity %f -> %d out of %d \n', sp_/length(wcoef), sp_, length(wcoef));

                try
                    figure((i-1)*3+1); clf; set(gcf, 'Position', [0 0 800 800]);
                    if ~use_real_visibilities
                        subplot(2,2,1); hold on; imagesc(log10(max(imn, 1e-8))); colorbar; axis image; title('Original image');
                    end
                    subplot(2,2,2); hold on; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
                    subplot(2,2,3); hold on; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');
                    subplot(2,2,4); hold on; imagesc(residualn_); colorbar, axis image; title('Residual image');
                end
                if ~isempty(aL2_v)
                    figure((i-1)*3+2); clf; set(gcf, 'Position', [0 0 800 400]);
                    if ~isempty(aL1_v)
                        hold on; plot(aL1_v/(min(aL1_v) + max(aL1_v)), 'r');
                    end
                    hold on; plot(aL2_v/(min(aL2_v) + max(aL2_v)), 'b');
                    hold on; plot(aepsilon/(min(aL2_v) + max(aL2_v)) * ones(length(aL2_v), 1), 'g');
                    xlabel('iterations');
                    ylabel('norm value (scaled)');
                    if ~isempty(aL1_v)
                        legend('L1 norm (scaled)', 'L2 norm (scaled)', 'L2 objective (scaled)');
                    else
                        legend('L2 norm (scaled)', 'L2 objective (scaled)');
                    end
                end

                if  ~isempty(aL1_vp) && ~isempty(aL2_vp)
                    figure((i-1)*3+3); clf; set(gcf, 'Position', [0 0 800 800]);
                    color = lines(R);
                    for q = 1:R
                        subplot(2,1,1); hold on; plot(aL2_vp(:, q), 'Color', color(q, :));
                        subplot(2,1,1); hold on; plot(aepsilonT{q} * ones(length(aL2_v), 1), 'Color', color(q, :), 'LineStyle', '--');
                    end
                    title('L2 norm for every node');
                    xlabel('iterations');
                    ylabel('norm value (scaled)');

                    color = lines(length(Psit));
                    for q = 1:length(Psit)
                        subplot(2,1,2); hold on; plot(aL1_vp(:, q), 'Color', color(q, :));
                    end
                    title('L1 norm for every node');
                    xlabel('iterations');
                    ylabel('norm value (scaled)');
                end
            end
        end

        if gen_only_average_figures
            
            if ~use_real_visibilties
                snr = snr/num_tests;
            end
            
            sp = round(sp/num_tests);
            soln = soln/num_tests;
            dirtyn = dirtyn/num_tests;
            residualn = residualn/num_tests;

            if ~use_real_visibilties
                fprintf('\n\n   Averege SNR: %f \n', snr);
            end
            fprintf('   Average solution sparsity %f -> %d out of %d \n', sp/length(wcoef), sp, length(wcoef));


            figure(1); clf; set(gcf, 'Position', [0 0 800 800]);
            if ~use_real_visibilties
                subplot(2,2,1); hold on; imagesc(log10(max(imn, 1e-8))); colorbar; axis image; title('Original image');
            end
            subplot(2,2,2); hold on; imagesc(log10(max(soln, 1e-8))), colorbar, axis image; title('Recovered image');
            subplot(2,2,3); hold on; imagesc(log10(max(dirtyn, 1e-8))), colorbar, axis image; title('Dirty image');
            subplot(2,2,4); hold on; imagesc(residualn), colorbar, axis image; title('Residual image');

        end

        clear wcoef;
        clear ys;
    end
