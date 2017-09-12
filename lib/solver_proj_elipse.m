function [x, k] = solver_proj_elipse(a, U, r, p, y, max_itr, min_itr, eps)
    r = r^2 - y' * y;
    U = [U; U].^2;
        
    a = [real(a); imag(a)];
    y = [real(y); imag(y)];
    p = [real(p); imag(p)];
    ysrt = sqrt(U) .* y;
    fa = a' * (U .* a) - 2 * a' * ysrt;

    ga = U .* a - ysrt;
    k = 0;
    if ( fa <= r )
        x = a;
    else
        x = p;
        u = U .* x - ysrt;
        v = a - x;
        while k < min_itr || (k < max_itr && (norm(u) == 0 || 1 - u' * v / (norm(u) * norm(v)) > eps))
            [x] = algo_4(k, 1, 1, 0.1, 0.8, x, a, U, fa, ga, eps, r, y, ysrt, u, v);
            u = U .* x - ysrt;
            v = a - x;
            k = k + 1;
        end
    end

    x = x(1:end/2) + 1i * x(end/2+1:end);
end

function [gamma, u, v] = algo_1(x, a, U, ~, ~, ~, ~, ~, ysrt)

    % algo 1
    v = a - x;
    vv = v' * v;

    u = U .* x - ysrt;
    z = u - u' * v / vv * v;

    zUz = z' * (U .* z);
    zz = z' * z;

    vUv = v' * (U .* v);
    vUz = v' * (U .* z);

    rho = 0.5 * (vUv / vv + zUz / zz + sqrt((vUv/vv - zUz/zz)^2 + 4 * vUz^2 / (vv * zz)));
    gamma = 1/rho;

end

function [gamma, u, v] = algo_2(~, ~, U, ~, ~, ~, ~, ~, ~, u, v)

    %% algo 2
            
%     u = U .* x - ysrt;
%     v = a - x;
    vv = v' * v;
    vUv = v' * (U .* v);

    z = u - u' * v / vv * v;
    zz = z' * z;
    zUz = z' * (U .* z);

    vUz = v' * (U .* z);
    zUv = z' * (U .* v);

    al = [sqrt(vv); 0];

    b = [u' * v / sqrt(vv); u' * z / sqrt(zz)];
    A = [vUv / vv vUz / (sqrt(vv)*sqrt(zz));
         zUv / (sqrt(vv)*sqrt(zz)) zUz / zz];

    [V_,S,V] = svd(A);
    s = diag(S);

    ah = V * al + 1./s .* (V * b);

    beta = (V * b)' * (1./s .* (V * b));
    l1 = s(1);
    l2 = s(2);

    pol = [l1 * (l1 - l2)^2; ...
            2 * l1 * l2 * ah(1) * (l1 - l2); ...
            l1^2 * l2 * ah(2)^2 + l1 * (l2 * ah(1))^2 - beta * (l1 - l2)^2; ...
            - 2 * beta * l2 * ah(1) * (l1 - l2); ...
            - beta * (l2 * ah(1))^2];
    rt = roots(pol);

    ahs = zeros(2, 1);

    for k = 1:4
        if imag(rt(k)) < 1e-15
            if (ah(1) - rt(k)) / (l1 * rt(k)) > 0
                ahs(1) = rt(k);
                ahs(2) = l1 * ah(2) * ahs(1) / ((l1 - l2) * ahs(1) + l2 * ah(1));
                break;
            end
        end
    end

    als = V_ * (ahs - 1./s .* (V * b));

    nu = 1 - als(1) / sqrt(vv) + u' * v / vv * als(2)/sqrt(zz);
    nu = real(nu);
    gamma = - als(2)/sqrt(zz) / nu;

end

function [x, u, v] = algo_3(k, x, a, U, fa, ga, eps, r, y)
    if mod(k, 3) == 0 || mod(k, 3) == 1
        [gamma, u, v] = algo_1(x, a, U, fa, ga, eps, r, y);

        c = x - gamma * u;
        w = c - a;
        wUw = w' * (U .* w);
        nu = - ga' * w / wUw - sqrt( (ga' * w / wUw)^2 - (fa - r) / wUw);
        nu = real(nu);
        x = a + nu * w;
    end

    if mod(k, 3) == 2
        [gamma, u, v] = algo_2(x, a, U, fa, ga, eps, r, y);

        c = x - gamma * u;
        w = c - a;
        wUw = w' * (U .* w);
        nu = - ga' * w / wUw - sqrt( (ga' * w / wUw)^2 - (fa - r) / wUw);
        nu = real(nu);
        x = a + nu * w;
    end
end


function [x, u, v] = algo_4(k, m1, m2, c1, c2, x, a, U, fa, ga, eps, r, y, ysrt, u, v)
     
    [gamma_2, u, v] = algo_2(x, a, U, fa, ga, eps, r, y, ysrt, u, v);

    if mod(k, m1+m2) < m1
        gamma = gamma_2;
        c = x - gamma * u;
        w = c - a;
        wUw = w' * (U .* w);
        nu = - ga' * w / wUw - sqrt( (ga' * w / wUw)^2 - (fa - r) / wUw);
        %nu = real(nu);
        x = a + nu * w;
    end

    if mod(k, m1+m2) >= m1
        [gamma_1, u, v] = algo_1(x, a, U, fa, ga, eps, r, y, ysrt);

        gamma = c1 * gamma_1 + c2 * gamma_2;
        c = x - gamma * u;
        w = c - a;
        wUw = w' * (U .* w);
        nu = - ga' * w / wUw - sqrt( (ga' * w / wUw)^2 - (fa - r) / wUw);
        %nu = real(nu);
        x = a + nu * w;
    end
end