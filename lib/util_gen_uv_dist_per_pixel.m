function [aY, aUV, aPos] = util_gen_uv_dist_per_pixel(u, v, Nox, Noy, y)

um = u;
vm = v;

aY = cell(Noy, Nox);
aUV = cell(Noy, Nox);
aPos = cell(Noy, Nox);


lsv = linspace(-pi, pi, Noy+1);
lsu = linspace(-pi, pi, Nox+1);
[v_, sv] = sort(vm);


for k = 1:Noy
    [sfv_l, sfv_h] = util_sort_find(v_, lsv(k), lsv(k+1), k<Noy, k==Noy);
    sfv = sv(sfv_l:sfv_h);
    if ~isempty(sfv) 
        [u_, su] = sort(um(sfv));
        y_ = y(sfv);
        for j = 1:Nox
            [sfu_l, sfu_h] = util_sort_find(u_, lsu(j), lsu(j+1), j<Nox, j==Nox);
            sfu = su(sfu_l:sfu_h);
            if ~isempty(sfu)
                aY{k, j} = y_(sfu);
                aUV{k, j} = [um(sfv(sfu)), vm(sfv(sfu))];
                aPos{k, j} = sfv(sfu);
            end
        end
    end
end

end


