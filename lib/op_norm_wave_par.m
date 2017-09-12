function val = op_norm_wave_par(w, A, At, im_size, tol, max_iter, verbose, numWorkers)
% computes the maximum eigen value of the compund 
% operator AtA
x = randn(im_size);
x = x/norm(x(:));
init_val = 1;

R = length(A);
util_create_pool(numWorkers);

for k = 1:max_iter

    y = cell(R, 1);
    parfor jj = 1:R
        y{jj} = At{jj}(w{jj}.^2 .* A{jj}(x));
    end

    x(:) = 0;
    for jj = 1:R
        x = x + y{jj};
    end

    val = norm(x(:));
    rel_var = abs(val-init_val)/init_val;
    if (verbose > 1)
        fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < tol)
       break;
    end
    init_val = val;
    x = x/val;

end

if (verbose > 0)
    fprintf('Norm = %e \n\n', val);
end

end

