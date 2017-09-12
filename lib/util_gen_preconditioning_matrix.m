function [aW] = util_gen_preconditioning_matrix(u, v, param)

if ~isfield(param, 'gen_uniform_weight_matrix'), param.gen_uniform_weight_matrix = 0; end
if ~isfield(param, 'uniform_weight_sub_pixels'), param.uniform_weight_sub_pixels = 1; end

um = u;
vm = v;

aWw = ones(length(vm), 1);

if param.gen_uniform_weight_matrix == 1
    Noy = param.uniform_weight_sub_pixels*param.Noy;
    Nox = param.uniform_weight_sub_pixels*param.Nox;
    
    lsv = linspace(-pi, pi, Noy+1);
    lsu = linspace(-pi, pi, Nox+1);
    [v_, sv] = sort(vm);
    
    for k = 1:Noy
        [sfv_l, sfv_h] = util_sort_find(v_, lsv(k), lsv(k+1), k<Noy, k==Noy);
        sfv = sv(sfv_l:sfv_h); 
        if ~isempty(sfv) 
            [u_, su] = sort(um(sfv));
            for j = 1:Nox
                [sfu_l, sfu_h] = util_sort_find(u_, lsu(j), lsu(j+1), j<Nox, j==Nox);
                sfu = su(sfu_l:sfu_h);
                if ~isempty(sfu)
                    aWw(sfv(sfu)) = length(sfu);
                end
            end
        end
    end
end

aW = 1./aWw;

end


