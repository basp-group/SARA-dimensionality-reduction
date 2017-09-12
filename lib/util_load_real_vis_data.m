function [y, uw, vw, sigma_noise, nW] = util_load_real_vis_data(visibility_file_name, param)

load(visibility_file_name)

% sigma_noise = 1.253*mad(y_V,0); %MAD estimator
% sigma_noise = sqrt(2);
sigma_noise = 1;

u = uvw(:,1);
v = uvw(:,2);
bmaxProj = max(sqrt(u.^2+v.^2));

if param.use_shift
   pos = param.shift_position;
   shift = exp(1i *(pos(1)*u+ pos(2)*v));
   y = y_I .* shift;
else
   y = y_I;
end
y = y .* weights;

% Setting imaging params
dl = param.pixel_size;

vw = v * pi / (bmaxProj * dl);
uw = u * pi / (bmaxProj * dl);

nW = weights;
Nm = length(uw);




%% subsample
if param.use_undersamplig
    Nmn = round(param.p * param.image_size_Nx * param.image_size_Ny);
    
    if strcmp(param.fu_undersampling_type, 'uniform')    
        if Nmn > Nm
            error('Can''t undersample the UV coverage: Not enough points');
        end
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            y(r) = [];
            nW(r) = [];
            Nm = Nm - length(r);
        end
    end
    
    if strcmp(param.fu_undersampling_type, 'gaussian')
        Nmo = length(uw);
        
        p_max = mvnpdf([0 0],[0 0], [param.fu_g_sigma 0; 0 param.fu_g_sigma]);
        prob = p_max * rand(Nmo, 1);

        sel = (prob > mvnpdf([uw vw],[0 0], [param.fu_g_sigma 0; 0 param.fu_g_sigma]));
        uw(sel) = [];
        vw(sel) = [];
        y(sel) = [];
        nW(sel) = [];
        
        fprintf('Selected %d points out of %d using a Gaussian profile \n', length(uw), Nmo);
        
        Nm = length(uw);
        if Nmn > Nm
            error('Can''t undersample the UV coverage. \n Not enough points left after generating the Gaussian profile: only %d points left after resampling %d points', Nm, Nmn);
        end
        
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            y(r) = [];
            nW(r) = [];
            Nm = Nm - length(r);
        end
        fprintf('Keeping only %d frequency points out of %d ... \n', Nmn, Nmo);
    end
    
    if strcmp(param.fu_undersampling_type, 'general-gaussian')
        Nmo = length(uw);
        mggdpdf = @(x, M)  exp(-(sum((x * diag(1./M)) .* x, 2)).^param.fu_ggd_beta/(2*param.fu_ggd_rho^param.fu_ggd_beta));
        
        p_max = mggdpdf([0 0], [param.fu_g_sigma param.fu_g_sigma]);
        prob = p_max * rand(Nmo, 1);

        sel = (prob > mggdpdf([uw vw], [param.fu_g_sigma param.fu_g_sigma]));
        uw(sel) = [];
        vw(sel) = [];
        y(sel) = [];
        nW(sel) = [];
        
        fprintf('Selected %d points out of %d using a generalised Gaussian profile \n', length(uw), Nmo);
        
        Nm = length(uw);
        if Nmn > Nm
            error('Can''t undersample the UV coverage. \n Not enough points left after generating the Gaussian profile: only %d points left after resampling %d points', Nm, Nmn);
        end
        
        while (Nm > Nmn)
            r = randi(Nm, Nm - Nmn, 1);
            r = unique(r);
            uw(r) = [];
            vw(r) = [];
            y(r) = [];
            nW(r) = [];
            Nm = Nm - length(r);
        end
        fprintf('Keeping only %d frequency points out of %d ... \n', Nmn, Nmo);
    end
end
