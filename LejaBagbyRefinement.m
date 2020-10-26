function [Am, Bm, Rm, info] = LejaBagbyRefinement(F, Z, Xi, varargin)

% LEJA_BAGBY_REFINEMENT Rational approximant for nonlinear eigenvalue
% problems [Am,Bm,Rm,info] = LEJABAGBYREFINEMENT(F, Z, XI, OPTS) is a
% rational approximant based on Leja-Bagby procedure. This is an auxiliary
% function called by NEP2RAT, it should not be called on its own. 

%This function accepts the following optional parameters:
%
%    opts:
%       - .dmax: the maximum degree of the approximant. Default is 20.
%       - .tol1: the tolerance of the first phase. Default is 1e-11.
%       - .tol2: the tolerance of the second phas. Default is 1e-13.
%       - .phase2: the second phase, which could be an exact search
%          ('exactLB'), a relaxed ('relaxedLB')search or nothing ('').
%          Default is 'exactLB'.
%       - .verbose: set the verbosity of the code (0,1,2). Default is 1.
%       - .cyclic: the choice of the poles. If 0, it uses the Leja-Bagby
%          method. If 1, it repeats the poles cyclically. Default is 0.
%
% It also returns additional information info

if nargin == 4
    opts = varargin{1};
else
    opts = struct();
end
Params = iParseInputs(F, Z, opts);

% Setting the parameters
dmax = Params.dmax;
tol1 = Params.tol1;
tol2 = Params.tol2;
Z2 = Params.Z2;
Z = union(Z, Z2); % we want Z being superset of Z2;
exactSearch = Params.exactSearch;
cyclic = Params.cyclic;
sparseFlag = Params.sparseFlag;
verbose = Params.verbose;
nF = Params.nF;
% End of setting the parameters

if any(isnan(Xi))
    error('All poles in Xi must be finite or infinite numbers.');
end

% If we plan to do an exact search, we pre compute all the F(Z(j))
maxF = nF;
FZ2 = cell(1, length(Z2));
if exactSearch
    FZ2 = arrayfun(F,Z2,'UniformOutput',false);
    % This is a constant to be used in second stopping criterion
    maxF = max(cellfun(@(z) norm(z, 'fro'), FZ2));
    if issparze(FZ2{1})
        % We save FZ2 in a matrix in the sparse case for a fast
        % implementation
       FZ2 = iDensify(FZ2);
    end
end

beta = [];
z = Z(1);
D{1} = F(z(1));
pol= Xi(1);
bZ = 1;
bXi = 1;
n = length(D{1});

startSearch = 0;
errvecMat = [];
if verbose >= 2
   fprintf('Starting the approximation.\n')
end
for j = 1:dmax
     if abs(pol(j)) < 1e-6
        warning('LEJA_BAGBY_REFINEMENT: One of the poles in the linearisation is very close to zero. Consider shifting the NLEP for stability.');
     end

    % get scaling parameter beta and next sigma
    bZ = bZ.*(Z - z(j))./(1 - Z/pol(j)); %b_j(Z)*beta_j
    [beta(j),ind] = max(abs(bZ)); % beta_j
    bZ = bZ/beta(j);
    bXi = bXi.*(Xi - z(j))./(1 - Xi/pol(j))/beta(j); %b_j(Xi)
    if startSearch
        % If we are in the second phase, we look for the next support point
        % considering the matrix valued function. So we overwrite the ind
        % of the instruction     [beta(j),ind] = max(abs(bZ));
        [err, ind] = iSearchNextPoint(FZ2, Rd, Z,  tol2*maxF, exactSearch);
        errvecMat = [errvecMat; err];
        if ind == 0
            % we did not find a point where the matrix error is greater
            % than the tolerance, so we have finished. We removed the last
            % element of errvecMat because it is the default value of
            % iSearchNextPoint
            errvecMat = errvecMat(1:end-1);
            % Normalize errvecMat
            errvecMat = errvecMat/maxF;
            break
        end
    end
    z(j+1) = Z(ind); % z_j
     if cyclic % cyclic poles
        pol(j+1) = Xi(mod(j,length(Xi))+1);
    else % leja-bagby
        [~,ind] = min(abs(bXi));
        pol(j+1) = Xi(ind);
    end

    % Rational basis functions.
    b = @(j, Lambda) arrayfun(@(lambda) ...
    prod((lambda-z(1:j))./(beta(1:j).*(1-lambda./pol(1:j)))),Lambda);
    % Next divided difference.
    Fz = F(z(j+1));
    % if we don't have a second phase, we need to compute maxF
    maxF = max(maxF, norm(Fz, 'fro')/sqrt(n));
    for k = 0:j-1
        Fz = Fz - b(k,z(j+1))*D{k+1};
    end
    D{j+1} = Fz/b(j,z(j+1));

    if ~startSearch
            errLB(j) = norm(D{j+1},'fro')/maxF;
            kappa = 3;
        if errLB(j) < tol1/kappa
            if exactSearch
                if verbose >= 2
                    if exactSearch == 1
                        fprintf('Starting the exact search of the points.\n')
                    else
                        fprintf('Starting the relaxed search of the points.\n')
                    end
                end
                startSearch = 1;
                Rd = @(zz) iRmMatrixHandleLejaBagbyFast(zz, b, D);
            else
        % we don't include this last point: in this way, if the first
        % divided difference matrix already satisfy the threshold, we don't
        % add an additional point
                z(j+1) = [];
                pol(j+1) = [];
                beta(j) = [];
                D(j+1) = [];
          %      errLB(j) = [];
                Rd = @(zz) iRmMatrixHandleLejaBagbyFast(zz, b, D);
                break
            end
        end
    else % We are already searching with the matrix, so we update R(z)
        Rd = @(zz) iRmMatrixHandleLejaBagbyFast(zz, b, D);
    end

end

% We need to define it here as well in the corner case the algorithm stays
% in phase 1, but it does not reach the approximation error before the
% maximum amount of iterations. This also covers the case of going into the
% second phase in the sparse case.
Rd = @(zz) iRmMatrixHandleLejaBagby(zz, b, D);

% Build the phase vector
phase = [ones(length(errLB),1); 2*ones(length(errvecMat),1)];

if verbose >= 2
   fprintf('Starting the linearization.\n')
end

[Am, Bm] = iLinearizeLejaBagby(z, pol, beta, D, sparseFlag);
% Setting the info parameter
info.Z2 = Z2;
info.zLB = z(:);
info.pol = pol(:);
info.phase = phase;
info.degree = length(z)-1;
info.errvec = errLB(:);
info.errvec2 = errvecMat(:);
info.D = D;
% Return a warning messagge if we reach the maximum degree
if info.degree == dmax
   if isempty(find(phase == 2, 1))
    if verbose >= 1
        warning(['The algorithm reached the maximum number of iterations', ...
            ' during the scalar Leja-Bagby portion.'])
    end
   info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the scalar Leja-Bagby portion.'];
   else
      if verbose >= 1
           warning(['The algorithm reached the maximum number of iterations', ...
            ' during the matrix Leja-Bagby portion.'])
      end
      info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the matrix Leja-Bagby portion.'];
   end
end

Rm = Rd;

end


%% Auxiliary functions

function [Am, Bm] = iLinearizeLejaBagby(z, xi, beta, D, sparseFlag)

% A simple implementation of the NLEIGS linearization by Guttel, Van
% Beeumen, Meerbergen and Wim.

z = z(:);
xi = xi(:);
beta = beta(:);
n = size(D{1}, 1); % we assume matrix square
N = length(D)-1;

if N == 0
   % we return a warning
   warning('We do not return a 0th linearization.')
   Am = [];
   Bm = [];
   return
end

if sparseFlag
   I = speye(n,n);
else
    I = eye(n,n);
end
   % Build the diagonal and subdiagonal of Am
   aux = diag([1; beta(1:N-1)]) + diag(z(1:N-1), -1);
   Am = kron(aux, I);
   % Build first block row
   Am(1:n, :) = cell2mat(D(1:N));
   Am(1:n, end-n+1:end) = Am(1:n, end-n+1:end)-z(N)*D{N+1}/beta(N);
   % Build the diagonal and subdiagonal of Bm
   aux1 = [1; beta(1:N-1)./xi(1:N-1)];
   aux = diag(aux1) + diag(ones(N-1,1), -1);
   Bm = kron(aux, I);
   % Build first block row
   Bm(1:n, :) = cell2mat(D(1:N))/xi(N);
   Bm(1:n,end-n+1:end) = Bm(1:n,end-n+1:end) - D{N+1}/beta(N);

end

function Rd = iRmMatrixHandleLejaBagbyFast(lambda, b, D)
% Build a matrix function in the divided differences form. If the matrices
% in the cell D are sparse, then we return a compact (vector)
% representation of Rd, which is much faster to evaluate.
if issparse(D{1})
    ind = find(D{randperm(length(D),1)});
    Rd = zeros(length(ind), length(D));
    bvec = zeros(length(D),1);
    for k = 1:length(D)
    Rd(:,k) = D{k}(ind);
    bvec(k) = b(k-1, lambda);
    end
    Rd = Rd*bvec;
else
    Rd = b(0, lambda)*D{1};
    for k = 1:length(D)-1
        Rd = Rd + b(k, lambda)*D{k+1};
    end
end
end

function Rd = iRmMatrixHandleLejaBagby(lambda, b, D)
% Build a matrix function in the divided differences form.
Rd = b(0, lambda)*D{1};
for k = 1:length(D)-1
    Rd = Rd + b(k, lambda)*D{k+1};
end
end

function [maxErr, ind] = iSearchNextPoint(FZ, RRD, Z,  tol, exactSearch)
% This is an auxiliary function to do either an exact search or relaxed
% search of the next point in the AAA algorithm in the case of a matrix
% function.
l = length(Z);
ind = 0;
maxErr = tol;
% Usually maxAux and MaxErr will be the same. We cover the corner case
% where the error is < tol to return max error, instead of the tolerance.
maxAux = 0;
if iscell(FZ)
    %FZ is a cell, because the function F(z) is dense.
    if exactSearch == 1
        for j = 1:l
            aux = norm(FZ{j} - RRD(Z(j)), 'fro');
            maxAux = max(aux, maxAux);
            if aux > maxErr
                maxErr = aux;
                ind = j;
            end
        end
    else % exactSearch == 0.5
        % We do a randomized shuffle to avoid a bias on the initial points
        indeces = randperm(length(Z));
        for k = 1:l
            j = indeces(k);
            aux = norm(FZ{j} - RRD(Z(j)), 'fro');
            maxAux = max(aux, maxAux);
            if aux > maxErr
                maxErr = aux;
                ind = j;
                break
            end
        end

    end
else
    % FZ is a matrix, coming from the iDensify when F is a sparse function;
    % Similarly, RRD(Z(j)) is a vector
    if exactSearch == 1
        for j = 1:l
            aux = norm(FZ(:,j) - RRD(Z(j)), 'fro');
            maxAux = max(aux, maxAux);
            if aux > maxErr
                maxErr = aux;
                ind = j;
            end
        end
    else % exactSearch == 0.5
        % We do a randomized shuffle to avoid a bias on the initial points
        indeces = randperm(length(Z));
        for k = 1:l
            j = indeces(k);
            aux = norm(FZ(:,j) - RRD(Z(j)), 'fro');
            maxAux = max(aux, maxAux);
            if aux > maxErr
                maxErr = aux;
                ind = j;
                break
            end
        end

    end
end
if ind == 0
   % it means that maxErr is still = tol, therefore we update it
   maxErr = maxAux;
end

end

function FzDense = iDensify(FzCell)
% This auxiliary function takes as an input a cell of length n of sparse
% matrices and return a dense matrix of size ind-by-n, where ind is the
% number of nonzero eleemnts in each sparse matrix of the cell

ind = find(FzCell{randperm(length(FzCell),1)});

FzDense = zeros(length(ind), length(FzCell));
for j = 1:length(FzCell)
FzDense(:,j) = FzCell{j}(ind);
end

end


%% Parse inputs

function Params = iParseInputs(F, Z, opts)
% Settings default parameters
dfVal = iDefaultValues(Z);
Params.phase2 = dfVal.phase2;
Params.tol1 = dfVal.tol1;
Params.tol2 = dfVal.tol2;
Params.Z2 = dfVal.Z2;
Params.dmax = dfVal.dmax;
Params.cyclic = dfVal.cyclic;
Params.verbose = dfVal.verbose;
Params.nF = dfVal.nF;
% End of default


if isfield(opts, 'phase2')
    Params.phase2 = opts.phase2;
end
if isfield(opts, 'tol1')
    Params.tol1 = opts.tol1;
end
if isfield(opts, 'tol2')
    Params.tol2 = opts.tol2;
end
if isfield(opts, 'Z2')
    Params.Z2 = opts.Z2(:);
end
if isfield(opts, 'dmax')
    Params.dmax = opts.dmax;
end
if isfield(opts, 'cyclic')
   Params.cyclic = opts.cyclic;
end
if isfield(opts, 'verbose')
    Params.verbose = opts.verbose;
end
if isfield(opts, 'nF')
    Params.nF = opts.nF;
end

% Check sparsity
if issparse(F(1))
    Params.sparseFlag = 1;
else
    Params.sparseFlag = 0;
end

%Convert Params.phase2 in the correspective easier parameters to parse
if isequal(Params.phase2,'exactLB')
    Params.exactSearch = 1;
elseif isequal(Params.phase2, 'relaxedLB')
    Params.exactSearch = 0.5;
elseif isequal(Params.phase2, '') || isequal(Params.phase2, 'LB')
    Params.exactSearch = 0;
else
    error(['Input in opts.phase2 not recognized.'])
end



end

function DefaultValues = iDefaultValues(Z)

DefaultValues.phase2 = 'exact';
DefaultValues.tol1 = 1e-11;
DefaultValues.tol2 = 1e-13;
DefaultValues.Z2 = Z;
DefaultValues.dmax = 20;
DefaultValues.cyclic = 0;
DefaultValues.verbose = 1;
DefaultValues.nF = 0;

end
