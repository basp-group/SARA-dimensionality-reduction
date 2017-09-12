
%%
% try
%     [1 2] * [1 1]
% 
%     %%
%     for k=1:1:100
% 
%         figure(k); subplot 121; imagesc(log10(abs((At(Gw'*(y{1}{1}-Gw*A(squeeze(result_st.sol_reweight_v{1}(k, :, :))))))))); axis image; colorbar
%         caxis([0.5 4]);
%         subplot 122; imagesc(log10(squeeze(result_st.sol_reweight_v{1}(k, :, :)))); axis image; colorbar
%         caxis([-5 0]); 
%     %     figure(k+200); 
%     %     histfit(col(real(At(Gw'*(y{1}{1}-Gw*A(squeeze(result_st.sol_reweight_v{1}(k, :, :))))))), 200); 
%         xcorr(k) = norm(corr2(abs(At(Gw'*(y{1}{1}-Gw*A(squeeze(result_st.sol_reweight_v{1}(k, :, :)))))), squeeze(result_st.sol_reweight_v{1}(k, :, :))));
%     end
% 
% end

%%


%%
% try

    %%
    R = 1
    for k=1:1:30

        figure(k); 
        set(gcf, 'Position', [50 1000 1200 300]);
        tt1 = real((At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :)))))));
        subplot 141; imagesc(tt1/norm(tt1(:))); axis image; colorbar
        caxis([-0.01 0.01]);
        subplot 142; imagesc(log10(squeeze(result_st.sol_best_bound_v{1}(k, :, :)))); axis image; colorbar
        caxis([-2.7 -1.5]); 
        subplot 143;
        im_res = col(real(At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :)))))));
        histfit(im_res, 200);
        title(sprintf('k = %f, s = %f', kurtosis(im_res), skewness(im_res)));
        subplot 144;
        vis_res = (cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))));
        vis_res = [real(vis_res); imag(vis_res)];
        histfit(vis_res, 200);
        title(sprintf('k = %f, s = %f', kurtosis(vis_res), skewness(vis_res)));
%         
%         No = length(W{1});
%         dirty = zeros(No, 1);
%         ns = A(squeeze(result_st.sol_best_bound_v{1}(k, :, :)));
%         for q = 1:R
%             dirty(W{q}) = dirty(W{q}) + G{q}' * ...
%                 (y{1}{q}-G{q}*ns(W{q}));
%         end
%         dirty = real(At(dirty));
        
%         tt1 = abs(real((At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))))))));
%         tt1 = abs(At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))))).^2);
        tt1 = At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :)))));
        tt2 = squeeze(result_st.sol_best_bound_v{1}(k, :, :));
        
%         tt1 = abs(tt1);
%         t1 = tt1(:);%'/norm(tt1(:));
%         t2 = tt2(:);
%         t1 = t1 - mean(t1);
%         t2 = t2 - mean(t2);
%         ixcr(k) = (t1(:)' * t2(:)) / norm(t1(:)) / norm(t2(:));
        
        tt1 = abs(real(tt1));
        ixcr(k) = norm(corr2(tt1/norm(tt1(:)), tt2)); %/norm(tt1(:))
    
        tx1 = real(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))));
        tx2 = real(Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))));
        vxcr(k) = norm(xcorr(tx1, tx2));
        
    end

% end




    %%
    R = 1
    for k=1:1:50

%         figure(k); 
%         set(gcf, 'Position', [50 1000 1200 300]);
%         tt1 = real((At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_v{1}(k, :, :)))))));
%         subplot 141; imagesc(tt1/norm(tt1(:))); axis image; colorbar
%         caxis([-0.01 0.01]);
%         subplot 142; imagesc(log10(squeeze(result_st.sol_v{1}(k, :, :)))); axis image; colorbar
%         caxis([-2.7 -1.5]); 
%         subplot 143;
%         im_res = col(real(At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_v{1}(k, :, :)))))));
%         histfit(im_res, 200);
%         title(sprintf('k = %f, s = %f', kurtosis(im_res), skewness(im_res)));
%         subplot 144;
%         vis_res = (cell2mat(y{1})-Gw*A(squeeze(result_st.sol_v{1}(k, :, :))));
%         vis_res = [real(vis_res); imag(vis_res)];
%         histfit(vis_res, 200);
%         title(sprintf('k = %f, s = %f', kurtosis(vis_res), skewness(vis_res)));
%         
%         No = length(W{1});
%         dirty = zeros(No, 1);
%         ns = A(squeeze(result_st.sol_best_bound_v{1}(k, :, :)));
%         for q = 1:R
%             dirty(W{q}) = dirty(W{q}) + G{q}' * ...
%                 (y{1}{q}-G{q}*ns(W{q}));
%         end
%         dirty = real(At(dirty));
        
%         tt1 = abs(real((At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))))))));
%         tt1 = abs(At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_best_bound_v{1}(k, :, :))))).^2);
        tt1 = At(Gw'*(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_v{1}(k, :, :)))));
        tt2 = squeeze(result_st.sol_v{1}(k, :, :));
        
%         tt1 = abs(tt1);
%         t1 = tt1(:);%'/norm(tt1(:));
%         t2 = tt2(:);
%         t1 = t1 - mean(t1);
%         t2 = t2 - mean(t2);
%         ixcr(k) = (t1(:)' * t2(:)) / norm(t1(:)) / norm(t2(:));
        
        tt1 = abs(real(tt1));
        ixcr(k) = sum(corr2(tt1/norm(tt1(:)), tt2)); %/norm(tt1(:))
    
%         tx1 = real(cell2mat(y{1})-Gw*A(squeeze(result_st.sol_v{1}(k, :, :))));
%         tx2 = real(Gw*A(squeeze(result_st.sol_v{1}(k, :, :))));
%         vxcr(k) = norm(xcorr(tx1, tx2));
        
    end

% end


%%
% y0f__ = A(result_st.sol{1});
% y0__ = cell(R, 1);
% for q = 1:R
%     y0__{q} = T{q} * y0f__(W{q});
% end

