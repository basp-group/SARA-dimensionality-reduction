function val = pow_method_par(A, At, im_size, R, tol, max_iter, verbose)
%Computes the maximum eigen value of the compund 
%operator AtA
%   
x=randn(im_size);
x=x/norm(x(:));
init_val=1;

for k=1:max_iter
    
    x1=At{1}(A{1}(x));
    for h = 2:R
        x1 = x1 + At{h}(A{h}(x));
    end
    val=norm(x1(:));
    rel_var=abs(val-init_val)/init_val;
    if (verbose > 0)
        fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < tol)
        break;
    end
    init_val=val;
    x=x1/val;
    
end


end

