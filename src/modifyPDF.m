function [new_pdf, alpha] = modifypdf(pdf, nb_meas)

% modifypdf checks PDF of the sampling profile and normalize it.
% Inputs
% pdf : sampling profile
% nb_meas : number of measurements.

if sum(pdf(:)<0)>0
    error('PDF contains negative values');
end
pdf = pdf/max(pdf(:));

% Find alpha
alpha_min = -1; alpha_max = 1;
while 1
    alpha = (alpha_min + alpha_max)/2;
    new_pdf = pdf + alpha;
    new_pdf(new_pdf>1) = 1; new_pdf(new_pdf<0) = 0;
    M = round(sum(new_pdf(:)));
    if M > nb_meas
        alpha_max = alpha;
    elseif M < nb_meas
        alpha_min = alpha;
    else
        break;
    end
end

