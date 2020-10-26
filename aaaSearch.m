function [r,pol,res,zer,z,f,w,errvec, errvecMat] = aaaSearch(FUN,F,Z,FZ2,Z2,tol1, tol2, mmax, exactSearch, verbose)
% Modification of aaa rational approximation of data F from a matrix
% function FUN on a set Z. [r,pol,res,zer,z,f,w,errvec, errvecMat] =
% aaaSearch(FUN,F,Z,tol,mmax, exactSearch). This is an AUXILIARY FUNCTION.
% It should NOT BE CALLED on its own. See HELP NEP2RAT for more
% information.


% Input: FUN: matrix function handle.
% F = vector of scalar data values, coming from u'*FUN(z)*v, where u and v
% are random vectors.
% Z = vector of sample points
% tol = relative tolerance tol, set to 1e-13 if omitted
% mmax: max type is (mmax-1,mmax-1), set to 100 if omitted
% exactSearch: flag for exact or relaxed search. See HELP HYBRID_SOLVER for
% more.
%
% Output: r = AAA approximant to F (function handle)
% pol,res,zer = vectors of poles, residues, zeros
% z,f,w = vectors of support pts, function values, weights
% errvec = vector of errors at each step
% errvecMat: empty vector if exactSearch is 0; otherwise, error of the
% matrix approximation at each iteration of the refinement


M = length(Z); % number of sample points
if nargin < 6, tol1 = 1e-13; end % default relative tol 1e-13
if nargin < 7, tol2 = 1e-15; end % default relative tol 1e-13
if nargin < 8, mmax = 100; end % default max type (99,99)
if nargin < 9, exactSearch = 0; end % default: no exactSearch
if ~isfloat(F), F = F(Z); end % convert function handle to vector
Z = Z(:); F = F(:); % work with column vectors
SF = spdiags(F,0,M,M); % left scaling matrix
J = 1:M; z = []; f = []; C = []; % initializations
errvec = []; R = mean(F);

%% % max(norm(F(Z2(j)), 'fro')), used as a relative tolerance ||F||_Z2
%% normFFZ = max(cellfun(@(zz) norm(zz, 'fro'), FZ2));

FzDense = []; % Tokens for the fast implementation
indDense = [];
sparseFlag = 0;
if issparse(FZ2{1})
   sparseFlag = 1;
   indDense = find(FZ2{randperm(length(FZ2),1)});
   FZ2 = iDensify(FZ2);
   FzDense = zeros(length(indDense),1);
end

startSearch = 0;
errvecMat = [];

for m = 1:mmax % main loop
    if ~startSearch
        [~,j] = max(abs(F-R)); % select next support point
    else
        % Next support point considering the matrix valued function
        [err, k] = iSearchNextPoint(FZ2, RRD, Z2, tol2, exactSearch); 
        errvecMat = [errvecMat; err];
        if k == 0
        % we did not find a point where the matrix error is greater than
        % the tolerance, so we have finished. We removed the last element
        % of errvecMat because it is the default value of iSearchNextPoint
        %errvecMat = errvecMat(1:end-1);
           break
        end
        % k is the index of the next point in the Z2 vector. We want the
        % index in the Z vector, which we know being a superset of Z2, thus
        % the next line works
        j = find(Z==Z2(k));
    end
    z = [z; Z(j)]; f = [f; F(j)]; % update support points, data values
    Fz{m} = FUN(Z(j)); %Update F(z_i)s
    if sparseFlag
        FzDense(:,m) = Fz{m}(indDense); 
    end
    J(J==j) = []; % update index vector
    C = [C 1./(Z-Z(j))]; % next column of Cauchy matrix
    Sf = diag(f); % right scaling matrix
    A = SF*C - C*Sf; % Loewner matrix
    [~,~,V] = svd(A(J,:),0); % SVD
    w = V(:,m); % weight vector = min sing vector
    N = C*(w.*f); D = C*w; % numerator and denominator
    R = F;
    R(J) = N(J)./D(J); % rational approximation
    
    if ~startSearch
        err = norm(F-R,inf);
        errvec = [errvec; err]; % max error at sample points
        if err <= tol1*norm(F,inf)
            if exactSearch
                if verbose >= 2
                    if exactSearch == 1
                        fprintf('Starting the exact search of the points.\n')
                    else
                        fprintf('Starting the relaxed search of the points.\n')
                    end
                end
                startSearch = 1;
                RRD = @(lambda) iRmatrixHandleFast(lambda,z,Fz,...
                    w, FzDense, indDense);
            else
                break
            end
        end %
    else % We are already searching with the matrix, so we update R(z)
        RRD = @(lambda) iRmatrixHandleFast(lambda,z,Fz, w, FzDense, ...
            indDense);
    end
    
end
r = @(zz) feval(@rhandle,zz,z,f,w); % AAA approximant as function handle
[pol,res,zer] = prz(r,z,f,w); % poles, residues, and zeros
%if ~exactSearch
[r,pol,res,zer,z,f,w] = ...
    cleanup(r,pol,res,zer,z,f,w,Z,F); % remove Frois. doublets (optional)
% %end

end


%% Auxiliary function

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

function Rd = iRmatrixHandle(zz,z, Fz,w)
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
    ind = find(Fz{randperm(length(Fz),1)});
    RdDense = zeros(length(ind), d);
    for j = 1:d
        RdDense(:,j) = Fz{j}(ind);
    end
    RdDense = RdDense*wizi.'/den;
    [m,n] = size(Fz{1});
    Rd = sparse(m,n);
    Rd(ind) = RdDense;
else
    Rd = wizi(1)*Fz{1};
    for j = 2:length(wizi)
        Rd  = Rd+wizi(j)*Fz{j};
    end
    Rd = Rd/den;
end
% Check whether Rd is NaN. If that is the case, it means zz = z_i for some i,
% therefore we force the interpolation.
if ~isempty(find(isnan(Rd),1))
    index = find(zz == z, 1);
    Rd = Fz{index};
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

function [pol,res,zer] = prz(r,z,f,w) % compute poles, residues, zeros
m = length(w); B = eye(m+1); B(1,1) = 0;
E = [0 w.'; ones(m,1) diag(z)];
pol1 = sort(eig(E,B)); 
pol = pol1(~isinf(pol1)); % finite poles to compute resids
pol = pol1(1:end-2); % we remove the spurious infinite poles, which are the last two;
indeces = isinf(pol);
pol(indeces) = 1e20; % we substitute real infinite poles with large numbers to avoid barycentric difficulties
dz = 1e-5*exp(2i*pi*(1:4)/4);
res = r(bsxfun(@plus,pol,dz))*dz.'/4; % residues
E = [0 (w.*f).'; ones(m,1) diag(z)];
zer = eig(E,B); zer = zer(~isinf(zer)); % zeros
end

function r = rhandle(zz,z,f,w) % evaluate r at zz
zv = zz(:); % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.'); % Cauchy matrix
r = (CC*(w.*f))./(CC*w); % AAA approx as vector
ii = find(isnan(r)); % find values NaN = Inf/Inf if any
for j = 1:length(ii)
    r(ii(j)) = f(find(zv(ii(j))==z)); % force interpolation there
end
r = reshape(r,size(zz)); % AAA approx
end

function [r,pol,res,zer,z,f,w] = cleanup(r,pol,res,zer,z,f,w,Z,F)
m = length(z); M = length(Z);
ii = find(abs(res)<1e-13); % find negligible residues
ni = length(ii);
if ni == 0, return, end
fprintf('%d Froissart doublets\n',ni)
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    % debug
    %abs(z(jj)-pol(ii(j)))/abs(pol(ii(j)))
    % end debug
    z(jj) = []; f(jj) = []; % remove nearest support points
end
for j = 1:length(z)
    F(Z==z(j)) = []; 
    Z(Z==z(j)) = [];
end
m = m-length(ii);
SF = spdiags(F,0,M-m,M-m);
Sf = diag(f);
C = 1./bsxfun(@minus,Z,z.');
A = SF*C - C*Sf;
[~,~,V] = svd(A,0); w = V(:,m); % solve least-squares problem again
r = @(zz) feval(@rhandle,zz,z,f,w);
[pol,res,zer] = prz(r,z,f,w); % poles, residues, and zeros
end
