function [Psi, Psit] = op_p_sp_wlt_basis(basis, nlevel, Ny, Nx)
% Resturns the operator to the sparsity wavelet basis passed as argument
% Each basis is considered to be distributed to a different node
%
% in:
% basis{:}    - cell of strigs with the names of the wavelets to be used
%             eg: {'db1', 'db2', 'self'}
% nlevel      - decomposition level
% Ny, Nx      - image size
% r_set       - flag to reset the sizes used for reconstruction
%
% out:
% Psi{:}[@]      - function handle for direct operator
% Psit{:}[@]     - function handle for adjoint operator


%% sparsity operator definition
dwtmode('per');
% construct a sting to repesent the desired inline function

Psit = cell(length(basis), 1);
for i = 1:length(basis)
    f = '@(x) [';
    if strcmp(basis{i}, 'self')
        f = sprintf('%s x(:);', f);
    else
        f = sprintf('%s wavedec2(x, %d, ''%s'')'';', f, nlevel, basis{i});
    end
    f = sprintf('%s]/sqrt(%d)', f, length(basis));
    Psit{i} = eval(f);
end


% for Psi it is a bit more complicated, we need to do some extra
% precomputations
Psi = make_Psi(basis, nlevel, Ny, Nx);

end

function Psi = make_Psi(basis, nlevel, Ny, Nx)
    Psi = cell(length(basis), 1);
    
    
    % estimate the structure of the data used to performe the
    % reconstruction
    S = cell(length(basis), 1);
    ncoef = cell(length(basis), 1);

    for i = 1:length(basis)
        if ~strcmp(basis{i}, 'self')
            [Cb, Sb] = wavedec2(zeros(Ny, Nx), nlevel, basis{i});
            S{i} = Sb;
            ncoef{i} = length(Cb(:));
        end
    end

    % construct a sting to repesent the desired inline function
    
    for i = 1:length(basis)
        f = '@(x)(';
        if strcmp(basis{i}, 'self')
            f = sprintf('%s reshape(x(%d:%d), [Ny Nx])', f, 1, Ny*Nx);
        else
            f = sprintf('%s waverec2(x(%d:%d), S{%d}, ''%s'')', f, 1, ncoef{i}, i, basis{i});
        end
        f = sprintf('%s)/sqrt(%d)', f, length(basis));
        Psi{i} = eval(f);
    end
end




