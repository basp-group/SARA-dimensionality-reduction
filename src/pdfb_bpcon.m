function xsol = pdfb_bpcon(y, epsilon, A, At, Psi, Psit, param)
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
sc = @(z) z*min(epsilon/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values
clip = @(z,T) sign(z).*min(abs(z),T); %clipping operator


%Initializations.

%Initial solution
if isfield(param,'initsol')
    xsol = param.initsol;
else
    xsol = 1/param.nu2*real(At(y)); 
end

%Initial dual variables
%L1 
r1 = Psit(xsol);

%L2 ball
r2 = A(xsol);

%Initial dual variables
if isfield(param, 'initv1')
    v1 = param.initv1;
else
    v1 = zeros(size(r1));
    %v1 = r1;
end

%Initial dual variables
if isfield(param, 'initv2')
    v2 = param.initv2;
else
    v2 = zeros(size(r2));
    %v2 = r2;
end

%Initial primal gradient
g = Psi(v1) + At(v2);

%Step sizes computation

%Step size for the dual variables
sigma1 = 1.0;
sigma2 = 1.0;

%Step size primal 
tau = 1.1/(sigma1*param.nu1 + sigma2*param.nu2);

%Initial objective value
fval = sum(param.weights(:).*abs(r1(:))); 
flag = 0;


%Main loop. Sequential.
for t = 1:param.max_iter
    
    %Primal update
    xsol = hardt(xsol - tau*g);
    
    %Dual variables update
    
    %L1 function update
    dummy = Psit(xsol); 
    r1 = v1 + sigma1*(2*dummy - r1);
    v1 = clip(r1,sigma1);
    r1 = dummy;
    %Objective value
    prev_fval = fval;
    fval = sum(param.weights(:).*abs(dummy(:)));
    
    %L2 ball projection update
    res = A(xsol);
    r2 = v2 + sigma2*(2*res - r2);
    v2 = r2 - y - sc(r2 - y);
    r2 = res;
    %Norm of resudual
    res1 = norm(res(:)- y(:));
    
    %Relative change of objective function   
    rel_fval = abs(fval - prev_fval)/fval;
    
    %Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm = %e, rel_fval = %e\n', ...
            fval, rel_fval);
        fprintf(' epsilon = %e, residual = %e\n\n', epsilon, res1);
    end
    
    %Global stopping criteria
    if (rel_fval < param.rel_obj && res1 <= epsilon*1.001)
        flag = 1;
        break;
    end
    
    %Update primal gradient
    g = Psi(v1) + At(v2);
          
end


%Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Final residual = %e\n', res1);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf(' Objective function = %e\n', fval);
        fprintf(' Relative variation = %e\n', rel_fval);
        fprintf(' Final residual = %e\n', res1);
        fprintf(' epsilon = %e\n', epsilon);
    end
end

end

