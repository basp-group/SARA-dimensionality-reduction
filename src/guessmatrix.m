function [matrix] = guessmatrix(diagonly, operator, matrixrows, matrixcols)
% Guesses the matrix corresponding to a given operator
% by operating on different delta vectors

if diagonly
    maxnonzeros = min(matrixrows, matrixcols);
    operdiag = zeros(maxnonzeros, 1);
else
    matrix = zeros(matrixrows, matrixcols); % BIG BIG BIG
end
for i=1:matrixcols
    deltacol = sparse(i, 1, 1, matrixcols, 1, 1);
    currcol = operator(deltacol);  % SLOW SLOW SLOW
    if diagonly
        if i > maxnonzeros
            break
        end
        operdiag(i) = currcol(i);
    else
        matrix(:,i) = currcol;
    end
    clear deltacol
end
if diagonly
    matrix = sparse(1:maxnonzeros, 1:maxnonzeros, operdiag, matrixrows, matrixcols, maxnonzeros);
end
end

