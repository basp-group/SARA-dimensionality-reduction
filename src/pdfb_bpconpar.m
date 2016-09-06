function xsol = pdfb_bpconpar(y, epsilon, A, At, Psi, Psit, param, R)
%
% sol = pdfb_bpcon(y, epsilon, A, At, Psi, Psit, param) solves:
%
%   min ||Psit x||_1   s.t.  ||y-A x||_2 <= epsilon and x>=0
%
%
% y contains the measurements. A is the forward measurement operator and
% At the associated adjoint operator. Psit is a sparfying transform and Psi
% its adjoint. PARAM a Matlab structure containing the following fields:
%
%   General parameters:
% 
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the objective value (default:
%   1e-4)
%       The algorithm stops if
%           | ||x(t)||_1 - ||x(t-1)||_1 | / ||x(t)||_1 < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%   - param.weights: weights (default = 1) for a weighted L1-norm defined
%       as sum_i{weights_i.*abs(x_i)}
%           
%
% Author: Rafael Carrillo
% E-mail: rafael.carrillo@epfl.ch
% Date: Feb. 15, 2015
%

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'nu1'), param.nu1 = 1; end
if ~isfield(param, 'nu2'), param.nu2 = 1; end
if ~isfield(param, 'weights'), param.weights = 1; end


% Useful functions for the projection
sc = @(z,radius) z*min(radius/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values
clip = @(z,T) sign(z).*min(abs(z),T); %clipping operator


%Initializations.

%Initial solution
if isfield(param,'initsol')
    xsol = param.initsol;
else
    xsol = 1/param.nu2*real(At{1}(y{1})); 
end

%Initial dual variables
%L1 variables
r1 = Psit(xsol);

%Initial dual variables
if isfield(param, 'initv1')
    v1 = param.initv1;
else
    v1 = zeros(size(r1));
end

%L2 ball variables
for h = 1:R
    r2{h} = A{h}(xsol);
end

%Initial dual variables
if isfield(param, 'initv2')
    v2 = param.initv2;
else
    for h = 1:R
        %v2{h} = zeros(size(y{h}));
        v2{h} = r2{h};
    end
end

%Initial primal gradient
g = Psi(v1);
for h = 1:R
    u{h} = At{h}(v2{h});
    g = g + u{h};
end

%Step sizes computation

%Step size for the dual variables
sigma1 = 1.0;
sigma2 = 1.0;

%Step size primal 
tau = 1.0/(sigma1*param.nu1 + sigma2*param.nu2);

%Initial objective value
fval = sum(param.weights(:).*abs(r1(:))); 
flag = 0;

res1 = zeros(R,1);
epsilont = epsilon(1)^2;
for h = 2:R
    epsilont = epsilont + epsilon(h)^2;
end
epsilont = sqrt(epsilont);


%Main loop. Sequential.
for t = 1:param.max_iter
    
    %Primal update
    prev_xsol = xsol;
    xsol = hardt(xsol - tau*g);
    prev_xsol = 2*xsol - prev_xsol;
    
    %Dual variables update
    
    %L1 function update
    r1 = Psit(prev_xsol); 
    v1 = v1 + sigma1*r1;
    v1 = clip(v1,sigma1);
    %Objective value
    prev_fval = fval;
    fval = sum(param.weights(:).*abs(r1(:)));
    
    %L2 ball projection update
    %Parallel for all blocks
    for k = 1:R
        r2{k} = A{k}(prev_xsol);
        v2{k} = v2{k} + sigma2*r2{k} - y{k};
        v2{k} = v2{k} - sc(v2{k},epsilon(k));
        u{k} = At{k}(v2{k});
        %Norm of residual
        res1(k) = norm(r2{k}(:)- y{k}(:));
    end
    
    %Relative change of objective function   
    rel_fval = abs(fval - prev_fval)/fval;
    
    %Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm = %e, rel_fval = %e\n', ...
            fval, rel_fval);
        fprintf(' epsilon = %e, residual = %e\n\n', epsilont, norm(res1));
    end
    
    %Global stopping criteria
    if (rel_fval < param.rel_obj && prod(res1 <= epsilon))
        flag = 1;
        break;
    end
    
    %Update primal gradient
    g = Psi(v1);
    for h = 1:R
        g = g + u{h};
    end
          
end


%Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Final residual = %e\n', norm(res1));
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Relative variation = %e\n', rel_fval);
        fprintf(' Final residual = %e\n', norm(res1));
        fprintf(' epsilon = %e\n', epsilont);
    end
end

end

