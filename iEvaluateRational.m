function Rm = iEvaluateRational(r,coeffs,z,sparseFlag)
reval = r(z);
d = length(reval);
[m,n] = size(coeffs{1});
if sparseFlag
    % We want to know the nonzero elements of Rm in order to allocate the
    % correct quantity of memory
    nonZero = max(cellfun(@(z) length(nonzeros(z)), coeffs));
    Rm = sparse([], [], [], m, n, nonZero);
else
    Rm = zeros(size(coeffs{1}));
end
for kk = 1:d
    Rm = Rm + coeffs{kk}*reval(kk);
end
end
