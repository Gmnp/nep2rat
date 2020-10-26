function [r, pol, res, zer, z, f, w, errvec, errvecMat] = aaa_svSearch(FUN, F, varargin)
%AAA   Modification of set-valued aaa ( original code by Lietaert and
% Meerbergen) that supports an exact (or relaxed) search.
%   It is an AUXILIARY function and it SHOULD NOT BE USED on its own. It is
%   called by weighted_aaaRef.m.


% Parse inputs:
[F, Z, Z2, M, nF, tol1, tol2, mmax, exactSearch, verbose] = ...
    parseInputs(F, varargin{:});

% Left scaling matrix:
SF = spdiags(F, 0:-M:-M*(nF-1), M*nF, M);

% Initialize values
F = F(:);
R = mean(F);
z = zeros(mmax,1);
f = zeros(mmax,nF);
ind = zeros(mmax,nF);
H = zeros(mmax,mmax-1);
S = zeros(mmax,mmax-1);

Q = zeros(M*nF,0);
C = zeros(M,0);

FzDense = []; % Tokens for the fast implementation
indDense = [];
if exactSearch
    % We evaluate the matrix functions at all the points if and only if we are
    % doing the exact or relaxed search
    FZ2 = arrayfun(FUN,Z2,'UniformOutput',false);
    % max(norm(F(Z2(j)), 'fro')), used as a relative tolerance
    normFFZ = max(cellfun(@(zz) norm(zz, 'fro'), FZ2));
    if issparse(FZ2{1})
        % Find the indeces of the nonzero elements
        indDense = find(FZ2{1});
        % In the sparse case, we densify to make the function evaluation
        % faster.
        FZ2 = iDensify(FZ2);
        FzDense = zeros(length(indDense),1);
    end
end

errvecMat = [];
startSearch = 0;
simpleUpdate = 1;
for m = 1:mmax
    if ~startSearch
        [~,loc] = max(abs(F-R)); % select next support point
    else
        % next support point considering the amtrix valued function
        [err, loc2] = iSearchNextPoint(FZ2, RRD, Z2,  tol2*normFFZ, ...
            exactSearch);
        errvecMat = [errvecMat; err/normFFZ];
        if loc2 == 0
            % we did not find a point where the matrix error is greater than
            % the tolerance, so we have finished. We removed the last element
            % of errvecMat because it is the default value of iSearchNextPoint
            errvecMat = errvecMat(1:end-1);
            break
        end
        % loc2 is the index of the next point in the Z2 vector. We want the
        % index in the Z vector, which we know being a superset of Z2, thus
        % the next line works
        loc = find(Z==Z2(loc2));
    end
    loc = mod(loc,M);
    ind(m,:) =  loc + (M*(loc==0):M:(nF-1+(loc==0))*M);  % Get indices of the z_i
    z(m) = Z(ind(m,1));                           % Add interpolation point
    f(m,:) = F(ind(m,:));                         % Add function values
    Fz{m} = FUN(z(m));                            % Add Function Matrix valued
    if exactSearch && issparse(FZ2{1})
        %The fast implementation for the sparse case
        FzDense(:,m) = Fz{m}(indDense);
    end
    C(:,end+1) = 1./(Z - z(m));                   % Get column of the Cauchy matrix.
    C(ind(1:m,1),m) = 0;                          % Set the selected elements to 0
    
    v = C(:,m)*f(m,:);                            % Compute the next vector of the basis.
    v = SF*C(:,m)-v(:);
    
    % Update H and S to compensate for the removal of the rows. This cannot
    % be done if at the previous iteration the choleski factorization
    % returned an error due to the matrix being not pos def for numerical
    % errors. If this is the case (simpleUpdate=0), then we compute the
    % weights through a big SVD.
    if simpleUpdate
    q = Q(ind(m,:),1:m-1);
    q = q*S(1:m-1,1:m-1);
    [Si,cholFlag] = chol(eye(m-1,m-1)-q'*q); % fixing this, GMNP 20/12/19
    end
    if cholFlag
        if simpleUpdate == 1 && verbose > 1 
            warning(['The Choleski factorization cannot be performed ', ...
                 'because the matrix is not numerically posdef. The ', ...
                'weights are computed with an SVD.'])
        end
        simpleUpdate = 0;
        Sf =  spdiags(f, 0:-m:-m*(nF-1), m*nF, m);
        C1 = C;
        A1 = SF*C1 - kron(speye(nF),C1)*Sf;   % Loewner matrix.
        J = [1:nF*M]';
        J(ind(1:m,:)) = [];
        % Solve least-squares problem to obtain weights:
        [~, ~, V] = svd(A1(J,:), 0);
        w = V(:,m);
    else
        H(1:m-1,1:m-1) = Si*H(1:m-1,1:m-1);
        S(1:m-1,1:m-1) = S(1:m-1,1:m-1)/Si;
        S(m,m) = 1;
        Q(ind(1:m,:),:) = 0;
        
        nv = norm(v);
        H(1:m-1,m) = Q'*v;
        H(1:m-1,m) = S(1:m-1,1:m-1)'*H(1:m-1,m);
        HH = S(1:m-1,1:m-1)*H(1:m-1,m);
        v = v-(Q*HH);
        H(m,m) = norm(v);
        % Reorthogonalization is necessary for higher precision
        it = 0;
        while (it < 3) && (H(m,m) < 1/sqrt(2)*nv)
            h_new = S(1:m-1,1:m-1)'*(Q'*v);
            v = v - Q*(S(1:m-1,1:m-1)*h_new);
            H(1:m-1,m) = H(1:m-1,m) + h_new;
            nv = H(m,m);
            H(m,m) = norm(v);
            it = it+1;
        end
        v = v/H(m,m);
        
        % Add v
        Q(:,end+1) = v;
        
        % Solve small least squares problem with H
        [~,~,V] = svd(H(1:m,1:m));
        w = V(:,end);
    end
    % Get the rational approximation
    N = C*bsxfun(@times,w,f(1:m,:));       % Numerator
    D = C*bsxfun(@times,w,ones(m,nF));     % Denominator
    R = N(:)./D(:);
    R(ind(1:m,:)) = F(ind(1:m,:));
    
    if ~startSearch 
        % we are still searching considering only the scalar functions
        err1 = norm(F-R,inf)*nF;
        % new and sharper bound
        FmRmatrix = reshape(abs(F-R), [M,nF]);
        err = sum(max(FmRmatrix));
        %err = sum(max(reshape(abs(F-R), [M,nF])));
        errvec(m) = err; % max error at sample points
        errvec1(m) = err1;
        if err <= tol1
            if exactSearch
                if verbose >= 2
                    if exactSearch == 1
                        fprintf('Starting the exact search of the points.\n')
                    else
                        fprintf('Starting the relaxed search of the points.\n')
                    end
                end
                startSearch = 1;
                RRD = @(lambda) iRmatrixHandleFast(lambda,z,Fz,w,FzDense,indDense);
            else
                break
            end
        end 
    else 
        % Searching with the matrix function, so we need to update R(z)
        RRD = @(lambda) iRmatrixHandleFast(lambda,z,Fz,w,FzDense,indDense);
    end
end

f = f(1:m,:);
w = w(1:m);
z = z(1:m);


% Note: When M == 2, one weight is zero and r is constant.
% To obtain a good approximation, interpolate in both sample points.
if ( M == 2 )
    z = Z;
    f = F;
    w = [1; -1];       % Only pole at infinity.
    w = w/norm(w);   % Impose norm(w) = 1 for consistency.
    errvec(2) = 0;
end

% Remove support points with zero weight:
I = find(w == 0);
z(I) = [];
w(I) = [];
f(I,:) = [];

% Construct function handle:
r = @(zz) reval(zz, z, f, w);

% Compute poles, residues and zeros:
[pol, res, zer] = prz(r, z, f, w);

% Remove Froissart doublets:
[r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F, nF);




end % of AAA()

%% Auxiliary functions
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

function Rd = iRmatrixHandleFast(zz,z, Fz,w, FzDense, indDense)
% This function takes the parameters of the surrogate AAA algorithm and
% returns the rational function Rd evaluated at the point zz. z and w
% are the sample points and the weights, Fz are the matrices F(z_i)
% evaluated at the sample points and stored in a cell.

d = length(z);
CC = 1./bsxfun(@minus,zz,z.'); % Cauchy matrix
den = CC*w; % denominator
wizi = w.'.*CC; % vectors of product w_i/(z-z_i), for the numerator

% Smart implementation in the sparse case
if issparse(Fz{1})
    Rd = FzDense*wizi.'/den;
    
    % Check whether Rd is NaN. If that is the case, it means zz = z_i for some i,
    % therefore we force the interpolation.
    if ~isempty(find(isnan(Rd),1))
        index = find(zz == z, 1);
        Rd = FzDense(:,index);
    end
else
    Rd = wizi(1)*Fz{1};
    for j = 2:length(wizi)
        Rd  = Rd+wizi(j)*Fz{j};
    end
    Rd = Rd/den;
    
    % Check whether Rd is NaN. If that is the case, it means zz = z_i for some i,
    % therefore we force the interpolation.
    if ~isempty(find(isnan(Rd),1))
        index = find(zz == z, 1);
        Rd = Fz{index};
    end
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


%% Evaluate rational function in barycentric form.
% The next functions were originally coded by the authors of set-valued
% aaa, a couple of bugs needed to be corrected
function r = reval(zz, zj, fj, wj)
% Evaluate rational function in barycentric form.
l = length(zz);
zv = zz(:);                             % vectorize zz if necessary
CC = 1./bsxfun(@minus, zv, zj.');       % Cauchy matrix
r = (CC*bsxfun(@times,wj,fj))./(CC*wj);             % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv),:) = kron(ones(sum(isinf(zv)),1),sum(bsxfun(@times,wj,fj),1)./sum(wj));

% Deal with NaN:
ii = find(isnan(r));
ii = [mod(ii(:)-1,length(zz)),floor((ii(:)-1)/l)]+1; % off-by-1 error corrected by GMNP 4/05/20
for jj = 1:size(ii,1)
    if ( isnan(zv(ii(jj,1))) || ~any(zv(ii(jj,1)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj,1),ii(jj,2)) = fj(zv(ii(jj,1)) == zj, ii(jj,2));
    end
end


end % End of REVAL().


%% Compute poles, residues and zeros.
% Function by GMNP 14/06/20
function [pol, res, zer] = prz(r, zj, fj, wj)
% Compute poles, residues, and zeros of rational function in barycentric form.
m = length(wj);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

% Compute residues via discretized Cauchy integral:
dz = 1e-5*exp(2i*pi*(1:4)/4);
nF = size(fj,2);

poldz = bsxfun(@plus, pol, dz);
rpoldz = r(poldz); % this is a matrix of size prod(size(poldz))-by-nF
res = zeros(size(poldz,1),nF);
for j = 1:nF
   aux = reshape(rpoldz(:,j), size(poldz)); %r_j(poldz)
   res(:,j) = aux*dz.'/4;
end

% Compute zeros via generalized eigenvalue problem:
for it = 1:size(fj,2)
    E = [0 (wj.*fj(:,it)).'; ones(m, 1) diag(zj)];
    zer{it} = eig(E, B);
    % Remove zeros of numerator at infinity:
    zer{it} = zer{it}(~isinf(zer{it}));
end
end % End of PRZ().


%% Cleanup
% Function by GMNP 14/06/20
function [r, pol, res, zer, z, f, w] = cleanup(r, pol, res, zer, z, f, w, Z, F, nF)
% Remove spurious pole-zero pairs.

% Reshape F as a matrix
Fm = reshape(F, [length(Z),nF]);

% Find negligible residues:
maxRes = max(res, [], 2); %max residual over the nF functions
ii = find(abs(maxRes) < 1e-13 * norm(F, inf));
ni = length(ii);
if ( ni == 0 )
    % Nothing to do.
    return
elseif ( ni == 1 )
    fprintf('1 Froissart doublet.\n')
else
    fprintf('%d Froissart doublets.\n', ni)
end

% For each spurious pole find and remove closest support point:
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    
    % Remove support point(s):
    z(jj) = [];
    f(jj,:) = [];
end

% Remove support points z from sample set:
for jj = 1:length(z)
    Fm(Z == z(jj),:) = [];
    Z(Z == z(jj)) = [];
end
m = length(z);
M = length(Z);

% Build Loewner matrix:
SF = spdiags(Fm, 0:-M:-M*(nF-1), M*nF, M);

Sf =  spdiags(f, 0:-m:-m*(nF-1), m*nF, m);
C = 1./bsxfun(@minus, Z, z.');      % Cauchy matrix.
A = SF*C - kron(speye(nF),C)*Sf;   % Loewner matrix.

% Solve least-squares problem to obtain weights:
[~, ~, V] = svd(A, 0);
w = V(:,m);

% Build function handle and compute poles, residues and zeros:
r = @(zz) reval(zz, z, f, w);
[pol, res, zer] = prz(r, z, f, w);

end % End of CLEANUP().


%% Parse Inputs:

function [F, Z, Z2, M, nF, tol1, tol2, mmax, exactSearch, verbose] = ...
    parseInputs(F, varargin)
% Input parsing for AAA.


% Sample points:
if ( ~isempty(varargin) && isfloat(varargin{1}) )
    % Z is given.
    Z = varargin{1};
    varargin(1) = [];
end

% Set defaults for other parameters:
tol1 = 1e-13;        % Relative tolerance.
tol2 = 1e-15;        % Relative tolerance for refinement.
mmax = 100;         % Maximum number of terms.

% Check if parameters have been provided:
while ( ~isempty(varargin) )
    if ( strncmpi(varargin{1}, 'tol1', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol1 = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'tol2', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            tol2 = varargin{2};
        end
        varargin([1, 2]) = [];
    elseif ( strncmpi(varargin{1}, 'Z2', 2) )
        if ( isfloat(varargin{2}))
            Z2 = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif ( strncmpi(varargin{1}, 'mmax', 4) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            mmax = varargin{2};
        end
        varargin([1, 2]) = [];
        
    elseif (strncmpi(varargin{1}, 'exactSearch', 11) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            exactSearch = varargin{2};
        end
        varargin([1,2]) = [];
   
    elseif (strncmpi(varargin{1}, 'verbose', 7) )
        if ( isfloat(varargin{2}) && isequal(size(varargin{2}), [1, 1]) )
            verbose = varargin{2};
        end
        varargin([1,2]) = [];
   
    else
        error('Argument unknown.')
    end
end

if ( exist('Z', 'var') )
    
    % Work with column vector:
    Z = Z(:);
    M = length(Z);
    
    % Function values:
    if ( isa(F, 'function_handle') || isa(F, 'chebfun') )
        % Sample F on Z:
        F = F(Z);
    elseif ( isnumeric(F) )
        % Work with column vector and check that it has correct length.
        if ( size(F,1) ~= M )
            error('Inputs F and Z must have the same length.')
        end
    else
        error('Input for F not recognized.')
    end
    
else
    Z = [];
    M = length(Z);
end
nF = size(F,2);

end % End of PARSEINPUT().


