function [approxErr,normFZ] = computeApproxErr(F, Rm, Z, varargin)
% COMPUTEAPPROXERR(F, Rm, norm) compute the relative error of the
% approximation R(z)~F(z) on the set Z. More precisely, it computes
%
%  (max_j norm(F(z_j) - Rm(z_j), val))/(max_j norm(F(z_j), val)),
%
% with z_j in Z. COMPUTEAPPROXERR accepts a optional parameters: "val" that
% orders which norm to use. The default value is 'fro'. The fifth parameter
% is RmSparse, which is used for sparse matrices when the norm to compute
% the error is not the Frobenius one. The sixth parameter is normFZ, which
% is an estimation of ||F||_Z = max_{z\in Z} ||F(z)||. The seventh
% parameter is tolnormest, which is used with the 2-norm when computing the
% 2-norm is too expensive and we substitute that command with normest

if nargin < 4 || isempty(varargin{1})
whichNorm = 'fro';
else
    whichNorm = varargin{1};
end
if nargin < 5 || isempty(varargin{2})
RmSparse = [];
else
   RmSparse = varargin{2};
end

if nargin < 6 || isempty(varargin{3})
normFZ = 0;
den = 0;
else
 normFZ = varargin{3};
 den = normFZ;
end

if nargin < 7 || isempty(varargin{4})
     tolnormest = 1e-6;
else
     tolnormest = varargin{4};
end

if whichNorm ~= 2
    if issparse(F(Z(1))) && ~issparse(Rm(Z(1)))
        num = 0;
        indDense = find(F(Z(1)));
        for k=1:length(Z)
            Rmz1 = Rm(Z(k));                                ,
            nnzRm = nnz(Rmz1);
            [~, ~, Fz1] = find(F(Z(k)));
            nnzFz = nnz(Fz1);
            if nnzRm == nnzFz
                num = max(num, norm(Fz1 - Rmz1, whichNorm));
            elseif nnzRm > nnzFz
                display('Used correction when sparsity patterns are different')
                num = max(num, norm(F(Z(k)) - RmSparse(Z(k)), whichNorm));
            else
                error('Different sparsity pattern between Rm and Fzm, which cannot be solved')
            end
            if normFZ==0, den = max(den, norm(Fz1, whichNorm)); end
        end
        approxErr = num/den;
        if normFZ==0, normFZ = den; end
    else
        num = 0;
        for k=1:length(Z)
            Fz1 = F(Z(k));
            Rmz1 = Rm(Z(k));
            num = max(num, norm(Fz1 - Rmz1, whichNorm));
            if normFZ==0, den = max(den, norm(Fz1, whichNorm));end
        end
        approxErr = num/den;
        if normFZ==0, normFZ = den; end
    end
else % if norm 2 and sparse, we have to use normest
    if issparse(F(Z(1)))
        num = 0;
        %    den = 0;
        for k=1:length(Z)
            Fz1 = F(Z(k));
            Rmz1 = Rm(Z(k));
            num = max(num, normest(Fz1 - Rmz1, tolnormest));
            if normFZ==0, den = max(den, normest(Fz1, tolnormest));end
        end
        approxErr = num/den;
        if normFZ==0, normFZ = den; end
        
    else
        num = 0;
        for k=1:length(Z)
            Fz1 = F(Z(k));
            Rmz1 = Rm(Z(k));
            num = max(num, norm(Fz1 - Rmz1, 2));
            if normFZ==0, den = max(den, norm(Fz1, 2));end
        end
        approxErr = num/den;
        if normFZ==0, normFZ = den; end
    end
end


end




