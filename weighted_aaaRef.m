function [Am, Bm, Rm, info] = weighted_aaaRef(F, Z, varargin)

% [Am, Bm, Rm, info] = WEIGHTED_AAAREF(F,Z,opts) rationally approximates
% F(z) with a variation of the set-valued AAA, named weighted AAA. This an
% AUXILIARY function called by NEP2RAT. See HELP NEP2RAT for more
% information.
%
% INPUT:
%  - F is the matrix valued function given as a structure. F.coeffs
% contains the matrix coefficients, while F.fun contains the scalar
% functions.
%  - Z is the set of point where the user wants the nonlinear function to
%  be approximated.

% Start retrieving parameters
Z = Z(:);
coeffs = F.coeffs;
d = length(coeffs);
[m,n] = size(coeffs{1});
fun = F.fun;
% End retrieving parameters

% Start optional parameters
if nargin == 3
    opts = varargin{1};
else
    opts = struct();
end
Params = iParseInputs(F, Z, opts);

tol1 = Params.tol1;
tol2 = Params.tol2;
Z2 = Params.Z2;
dmax = Params.dmax;
verbose = Params.verbose;
exactSearch = Params.exactSearch;
exactLB = Params.exactLB;
LBref = Params.LBref;
sparseFlag = Params.sparseFlag;

% End optional parameters

% Build function handle
F1 = @(zz) iEvaluateRational(F.fun, F.coeffs, zz, sparseFlag);

% Normalizing the coefficients
cnorms = zeros(1,d);
coeffsNorm = coeffs;
for kk = 1:d
    cnorms(kk) = norm(coeffs{kk}, 'fro');
    coeffsNorm{kk} = coeffs{kk}/cnorms(kk);
end
funNorm = @(z) fun(z).*cnorms;

if verbose >= 2
    fprintf('Using weighted_aaa on the scalar functions.\n')
end

%This is a lower bound of ||F||_{\Sigma}, which is used in the stopping
%criterion.
state = rng();
rng('default');
u = randn(m,1);
u = u/norm(u);
rng(state);
if sparseFlag
    uj = sparse(m,d);
else
    uj = zeros(m,d);
end
for jj = 1:d
    uj(:,jj) = coeffs{jj}*u;
end
normFSigma = 0;
for jj = 1:length(Z)
    normFSigma = max(normFSigma, norm(uj*fun(Z(jj)).'));
end
% End computation of lower bound of ||F||_Sigma
tol_aaa = tol1*normFSigma;
[r, pol, res, zer, z, f, w, errvec, errvecMat] = aaa_svSearch(F1, ...,
    funNorm, Z, 'Z2', Z2, 'tol1', tol_aaa, 'tol2', tol2, 'mmax', dmax+1, ...,
    'exactSearch', exactSearch, 'verbose', verbose);

if isempty(pol)
    pol = inf;
end

% Setting the output "info"
info.r = @(z) r(z)./cnorms;
info.pol = pol;
info.res = res;
info.zer = zer;
info.z = z;
info.f = f;
info.w = w;
info.Z2 = Z2;
info.errvec = errvec/max(cnorms); % error normalized
info.errvec2 = errvecMat(:);
phase = [ones(length(info.errvec),1); 2*ones(length(errvecMat),1)];
info.phase = phase;
info.degree = length(z)-1;
info.nF = normFSigma;
info.msg = '';
% End of setting the output "info"

% We want to return a warning if we reached the max degree of
% approximation. Due to the cleaning up of the Froissart doublets, the
% degree could be lower than the max degree, even if we reached dmax.
% Therefore the criterion with the vector phase.
if length(info.phase) == dmax+1
    if isempty(find(phase == 2, 1))
        if verbose >= 1
            warning(['The algorithm reached the maximum number of iterations', ...
                ' during the weighted AAA portion'])
        end
        info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the weighted AAA portion'];
    else
        if verbose >=1
            warning(['The algorithm reached the maximum number of iterations', ...
                ' during the weighted AAA with exact/relaxed search portion'])
        end
        info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the weighted AAA with exact/relaxed search portion'];
    end
end



% Simple linearization if refinement is not used
if ~(exactSearch || LBref) || (info.degree == dmax)
    if verbose >= 2
        fprintf('Linearizing the rational approximation.\n')
    end
    [Am, Bm] = iAAALinearize(z, f, w, coeffsNorm, sparseFlag);
    Rm = @(z) iEvaluateRational(r, coeffsNorm, z, sparseFlag);
    info.Rm = Rm;
    return
end

d = info.degree;

% matrices F(z_i), to build the rational approximation in barycentric form.
% They are stored in a cell of d+1 m-by-n matrices
Fz = arrayfun(F1,z,'UniformOutput',false);

% If we are doing an exact search in the Leja-Bagby refinement as the
% second phase, then we compute the norms of all the points in Z2
if exactLB
    FZ2 = arrayfun(F1,Z2,'UniformOutput',false);
    maxF = max(cellfun(@(z) norm(z, 'fro'), FZ2));
else
    FZ2 = cell(size(Z)); % empty cell to be used as a token in iLejaBagbyrefinement
    maxF = max(cellfun(@(z) norm(z, 'fro'), Fz));
end

% Build the matrix-handle of the approximation in the fast way: if the
% matrix is sparse, then it is densified.
FzDense = []; %tokens for fast implementation
indDense = [];
if sparseFlag
    indDense = find(Fz{1});
    FzDense = iDensify(Fz);
end
Rd = @(lambda) iRmatrixHandleFast(lambda,z,Fz,w. FzDense, indDense);
info.Rm = Rd;

% If we have done an exact search it means we have finished, so we need to
% build Am, Bm and return. Notice that we don't need to check LBref = 0,
% because in that case we would have already returned.
if exactSearch
    Rm = @(lambda) iRmatrixHandle(lambda,z,Fz,w);
    % iLinearize takes care of the linearization, both in the simple
    % weighted_aaa case and when the linearization is mixed due to
    % Leja-Bagby refinement. It distingueshes between the two using d
    % and dmax: when they are equal, it means it did not enter the
    % Leja-Bagby refinement, this is why we are giving two "ds" in
    % input. In addition, it would need the vector of beta and of
    % poles xi with the mixed linearization, so here they are empty.
    if verbose >= 2
        fprintf('Linearizing the rational approximation.\n')
    end
    [Am, Bm] = iLinearize(z,w,Fz,[],[],m,n,d,d, sparseFlag);
    return
end

%% Leja-Bagby Refinement
if verbose >= 2
    fprintf('Starting Leja-Bagby refinement.\n')
end

% We do not start the Leja-Bagby refinement if Z and pol have points in
% common.
if ~isempty(intersect(Z,pol))
    if verbose >= 1
        warning(['The Leja-Bagby refinement does not work if Z and Xi '...
            'have points in common.']);
    end
    Am = [];
    Bm = [];
    info.msg = 'Z and Xi have points in common.';
    Rm = @(z) 0;
    return
end

[Rm, z, w, Fz, beta, xi, errLB, RmFast] = iLejaBagbyRefinement(F1, FZ2, Rd, Z, Z2, ...
    z, w, pol, Fz, dmax, maxF, tol2, sparseFlag, exactLB);

% This is a densified function-handle when Rm is sparse. It is used in the
% tests to compute the approximation error in a fast way.
info.Rm = RmFast;

if verbose >=2
    fprintf('Linearizing the rational approximation.\n')
end
[Am, Bm] = iLinearize(z,w,Fz,beta,xi,m,n,d,dmax, sparseFlag);

%Update the info output
info.ordPol = xi;
info.zLB = z;
info.degree = length(z)-1;
phase = [phase; 2*ones(length(errLB),1)];
info.phase = phase;
info.errvec2 = errLB(:);
% Warning if we have reached the maximum degree of approximation
if info.degree == dmax
    if verbose >=1
        warning(['The algorithm reached the maximum number of iterations', ...
            ' during the Leja-Bagby refinement'])
    end
    info.msg = ['The algorithm reached the maximum number of iterations', ...
        ' during the Leja-Bagby refinement'];
end


end

%% Auxiliary functions

function [LE, LF] = iAAALinearize(z, f, w, coeffs, sparseFlag)
n = size(coeffs{1},1);
[m,d] = size(f);
% If it is large and not sparse, we make it sparse nonetheless to avoid
% swapping
if (m+1)*n > 1000
    sparseFlag = 1;
end

if sparseFlag
    F = -speye(m) + spdiags(ones(m-1,1),-1, m,m);
    E = -spdiags(z,0,m,m);
    if m > 1
        % In this way we cover the case when m = 1, where the second
        % addendum is 0. The dense case does not need this distinction.
        E = E + spdiags(z(1:m-1),-1,m,m);
    end
    In = speye(n);
    b = spalloc(m,1,1); b(1) = 1;
    LE = sparse(n*(m+1), n*(m+1));
    LF = LE;
else
    F = -eye(m) + diag(ones(1,m-1),-1);
    E = -diag(z) + diag(z(1:m-1), -1);
    In = eye(n);
    b = zeros(m,1); b(1) = 1;
    LE = zeros(n*(m+1), n*(m+1));
    LF = LE;
end
F(1,1) = 0;
E(1,:) = w;


a = f.*w; % the columns of a are the vectors a_i: f_i(z_j)*w_j
% Now build Etilde and Ftilde
Et = [-b E];
Ft = [zeros(m,1) F];

aux = kron(a(:,1).', coeffs{1});
for j = 2:d
    aux = aux + kron(a(:,j).', coeffs{j});
end
LE(1:n,n+1:end) = aux;
LE(n+1:end,:) = kron(Et, In);
LF(n+1:end,:) = kron(Ft, In);

end

function [Am, Bm] = iLinearize(z, w, Fz , beta, xi, m,n, d,dmax, mysparse)
% This auxiliary function linearize the rational approximation. We need two
% distinctions. The first one is caused by sparsity, because in the full
% case, Fz is a 3D array containing the matrices  of the evaluations F(z_i)
% and the divided differences matrices D_i, while in the sparse case, Fz is
% a cell. The second one is on the degree. If d=dmax, it means that we did
% not enter the Leja-Bagby refinement, therefore the linearization is a bit
% simpler.

% Another way we did not enter the Leja-Bagby refinement is if
% beta is empty. In order to avoid worst switch cases in other part of the
% code, if beta is empty, then we set dmax equal to d, which will start the
% correct linearization.
if isempty(beta)
    dmax = d;
end

if d == 0
    % we return a warning
    warning('We do not return a 0th linearization.')
    Am = [];
    Bm = [];
    return
end

if mysparse
    if d == dmax
        N = length(z)-1; % in this case, N = d
        I = speye(m,n);
        % Build the diagonal and subdiagonal part of Bm
        aux = diag(ones(d-1,1), -1) - eye(d);
        Bm = kron(aux, I);
        % Build the first row of Bm
        FFw = cell(1,d+1);
        for j = 1:d+1
            FFw{j} = Fz{j}*w(j); %F(z_i)*w_i
        end
        firstBlkRow = cell2mat(FFw)/z(end); %F(z_i)*w_i/z_d
        Bm(1:m,:) = firstBlkRow(:,1:end-n);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + ...
            firstBlkRow(:,end-n+1:end);
        % Build the diagonal of Am
        aux = diag(z(1:end-1)) + diag(z(1:end-2), -1);
        Am = kron(aux,I);
        % Build the first row of Am. We already have F(z_i)*w_i
        Am(1:m,:) = cell2mat(FFw(1:end-1));
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + ...
            z(end-1)*FFw{end}/z(end);
    else
        % we entered the refinement cycle, so the linearization is different
        beta = beta(:);
        N = length(z)-1; % k+l in the notes.pdf
        I = speye(m,n);
        % Build the diagonal and subdiagonal part of Bm
        aux = diag([ones(d, 1); 1./beta(1:end-1)], -1);
        aux = aux - diag(ones(N,1));
        Bm = kron(aux,I);
        % Build the first row of Bm
        FFw = Fz.';
        
        for j = 1:d+1
            FFw{j} = Fz{j}*w(j); %first d+1 term: %F(z_i)*w_i
        end
        firstBlkRow = cell2mat(FFw)/xi(end);
        Bm(1:m,:) = firstBlkRow(:,1:end-n);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + ...
            firstBlkRow(:,end-n+1:end)/beta(end);
        % Build the diagonal of Am
        aux1 = z(d+1:end-2)./beta(1:end-1);
        aux = diag([z(1:d); aux1], -1);
        aux = aux - diag([z(1:d+1); xi(d+1:end-1)]);
        Am = kron(aux,I);
        % Build the first row of Am
        Am(1:m,:) = cell2mat(FFw(1:end-1));
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + ...
            Fz{end}*z(end-1)/(beta(end)*xi(end));
    end
else % full matrix
    Fz = cell2mat(reshape(Fz, [1,1,length(Fz)]));
    if d == dmax
        N = length(z)-1; % in this case, N = d
        Am = zeros(m*N , n*N);
        Bm = Am;
        I = eye(m,n);
        % Build the diagonal and subdiagonal part of Bm
        aux = diag(ones(d-1,1), -1) - eye(d);
        Bm = kron(aux, I);
        % Build the first row of Bm
        FFw = bsxfun(@times, Fz, reshape(w, [1,1,d+1])); %F(z_i)*w_i
        firstBlkRow = bsxfun(@rdivide, FFw, z(end)); %F(z_i)*w_i/z_d
        Bm(1:m,:) = reshape(firstBlkRow(:,:,1:end-1), [m,n*N]);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + firstBlkRow(:,:,end);
        % Build the diagonal of Am
        aux = diag(z(1:end-1)) + diag(z(1:end-2), -1);
        Am = kron(aux,I);
        % Build the first row of Am. We already have F(z_i)*w_i
        Am(1:m,:) = reshape(FFw(:,:,1:end-1), [m,n*N]);
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + ...
            z(end-1)*FFw(:,:,end)/z(end);
    else
        % we entered the refinement cycle, so the linearization is different
        beta = beta(:);
        N = length(z)-1; % k+l in the notes.pdf
        Am = zeros(m*N , n*N);
        Bm = Am;
        I = eye(m,n);
        % Build the diagonal and subdiagonal part of Bm
        aux = diag([ones(d, 1); 1./beta(1:end-1)], -1);
        aux = aux -diag(ones(N,1));
        Bm = kron(aux,I);
        % Build the first row of Bm
        FFw = Fz;
        FFw(:,:,1:d+1) = bsxfun(@times, Fz(:,:,1:d+1), ...
            reshape(w, [1,1,d+1])); %F(z_i)*w_i
        firstBlkRow = bsxfun(@rdivide, FFw, xi(end));
        Bm(1:m,:) = reshape(firstBlkRow(:,:,1:end-1), [m,n*N]);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + ...
            firstBlkRow(:,:,end)/beta(end);
        
        % Build the diagonal of Am
        aux1 = z(d+1:end-2)./beta(1:end-1);
        aux = diag([z(1:d); aux1], -1);
        aux = aux - diag([z(1:d+1); xi(d+1:end-1)]);
        Am = kron(aux,I);
        % Build the first row of Am
        FFw = Fz;
        FFw(:,:,1:d+1) = bsxfun(@times, Fz(:,:,1:d+1), ...
            reshape(w, [1,1,d+1])); %F(z_i)*w_i
        Am(1:m,:) = reshape(FFw(:,:,1:end-1), [m,n*N]);
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + ...
            Fz(:,:,end)*z(end-1)/(beta(end)*xi(end));
    end
end


end



function [Rm, z, w, Fz, beta, xi, errLB, RmFast] = iLejaBagbyRefinement...
    (F, FZ2, Rd, Z, Z2, z, w, pol, Fz, dmax, maxF, tol2, mysparse, exactLB)

l = length(Z);
l2 = length(Z2);
d = length(z)-1;
dPol = length(pol);

bPol = ones(dPol,1);
bZ = ones(l,1);
%xi = zeros(size(pol));
xi = pol(:);
betaTemp = zeros(d,1);
betaTemp(1) = 1;
errvecMat = [];
RdStart = Rd;
RdSparse = @(zz) iRmatrixHandle(zz,z,Fz, w);
RdStartSparse = RdSparse; % if the matrices are dense, then RdSparse = Rd;
indDense = []; % token for fast implementation
FzDense = [];  % token
DDense = [];
if mysparse
    indDense = find(Fz{randperm(length(Fz),1)});
    FzDense = iDensify(Fz);
    DDense = zeros(length(indDense), 1);
    if exactLB
        FZ2 = iDensify(FZ2);
    end
end

% First cycle: reordering the poles returned by AAA.
for j = 1:d
    [~,ind] = min(abs((pol - z(j)).*bPol));
    if j < dPol+1
        xi(j,1) = pol(ind);
    else
        xi(j,1) = xi(mod(j-1,dPol)+1);
    end
    bZ = bZ.*(Z - z(j))./(1-Z/xi(j));
    betaTemp(j,1) = max(abs(bZ));
    bZ = bZ/betaTemp(j);
    bPol = bPol.*(pol-z(j))./(betaTemp(j).*(1-pol./xi(j)));
end

%ZM = Z \ z, to avoid barycentric form difficulties. Same for Z2M
ZM = setdiff(Z, z);
[Z2M, indexes]= setdiff(Z2,z);
if exactLB
    %matrix function F evaluated on set Z2M
    if mysparse
        FZ2M = FZ2(:, indexes); % FZ2 is a matrix
    else
        FZ2M = FZ2(indexes);    % FZ2 is a column cell
    end
end
CC = 1./bsxfun(@minus,ZM, z.'); % Cauchy matrix
bZM = (ZM-z(d+1)).^(-1)./(CC*w); % b_d(z) evaluated on ZM
beta(1) = 1;
%function b_{d}(z).
bd = @(zz) bfun(d, zz, z, w);

% Rd of degree d is fully formed here. d+1 interpolation conds and d poles
% Second cycle: Leja-Bagby refinement
for j = d+1:dmax
    % new pole chosen cyclically. No off-by-one because it starts from 1 in
    % the note as well.
    xi(j,1) = xi(mod(j-1,dPol)+1); 
    if abs(xi(j,1)) < inf
        bZM = bZM.*(ZM-z(j))./(1-ZM/xi(j)); % b_{j}(Z^{(j)})
    else
        bZM = bZM.*(ZM-z(j))./(1-ZM/xi(j)); % b_{j}(Z^{(j)})
    end
    %Compute z_{j+1}, beta_{j} and normalize. Again, notice the
    %off-by-one.
    [beta(j-d), ind] = max(abs(bZM));
    bZM = bZM/beta(j-d);
    % If we are going for an exact or relaxed search on the matrix
    % function, we look for another point Z(ind)
    if exactLB
        [err, indAux] = iSearchNextPoint(FZ2M, Rd, Z2M,  tol2*maxF, exactLB);
        errvecMat(j-d) =  err;
        if indAux == 0
            % we did not find a point where the matrix error is greater
            % than the tolerance, so we have finished. We removed the last
            % element of errvecMat because it is the default value of
            % iSearchNextPoint
            errvecMat = errvecMat(1:end-1);
            % Normalize errvecMat
            errvecMat = errvecMat/maxF;
            %We have computed one pole and one beta too many, so we pop
            %them out, otherwise the linearization breaks.
            beta(end) = [];
            xi(end) = [];
            break
        else
            % Remove the evaluation from FZ2M.
            if mysparse
                FZ2M(:,indAux) = [];  % FZ2M is a matrix
            else
                FZ2M(indAux) = [];  % FZ2M is a column cell
            end
            % indAux is the index of the new point with respect of the
            % vector Z2 and FZ2M. We want that same point, but with
            % the index of Z. Given that ZM is a superset of Z2M, the
            % next line works
            ind = find(ZM == Z2M(indAux));
            % Remove the point from Z2
            Z2M(indAux) = [];
        end
    end
    
    z(j+1,1) = ZM(ind);
    Fsigma = F(z(j+1));
    n = length(Fsigma); %should be passed on before
    D{j-d} = (Fsigma - RdSparse(z(j+1)) )/bZM(ind);
    % Stopping criterion suggested by Fran.
    nFF = norm(Fsigma, 'fro')/sqrt(n);
    if maxF<nFF, display('tighter lower bnd');end
    maxF = max(maxF, nFF);
    errLB(j-d) = norm(D{j-d}, 'fro')/maxF;
    kappa = 3;
    if errLB(j-d) < tol2/kappa
        % we don't include this last point: in this way, if the first
        % divided difference matrix already satisfy the threshold, we don't
        % add an additional point
        z(j+1) = [];
        xi(j) = [];
        beta(j-d) = [];
        D(j-d) = [];
      %  errLB(j-d) = [];
        break
    else
        Fz{j+1} = D{j-d};
        if mysparse
            FzDense(:,j+1) = Fz{j+1}(indDense);
            DDense(:,j-d) = D{j-d}(indDense);
        end
        ZM(ind) = [];
        % Remove the new point from b_{j+1}(Z^(j))
        bZM(ind) = [];
        % Update the approximation. Here we lose the ability to evaluate Rd on
        % arrays. First, we build b_{j+1}(z) from b_{j}(z). Again, z and beta are
        % off-by-one.
        b = @(k, Lambda) bd(Lambda)*arrayfun(@(lambda) ...
            prod((lambda-z(d+1:d+k))./(beta(1:k).'.*(1-lambda./xi(d+1:d+k)))),Lambda);
        Rd = @(zz) iRmatrixHandleLBFast(zz, RdStart, b, D, DDense);
        RdSparse = @(zz) iRmatrixHandleLB(zz, RdStartSparse, b, D);
    end
    
end

% Rational approximation. If the approximation is already good, it means we
% leave the cycle before computing any divided difference matrix, therefore
% we must return the rational approximation from the input. We can't simply
% return RdStart because it may be in the "densified" state if the function
% handle is sparse.

if (j == d+1)
    Rm = @(zz) iRmatrixHandle(zz,z,Fz, w);% Rm = RdStart
else
    Rm = @(zz) iRmatrixHandleLB(zz, RdStartSparse, b, D);
end

% Fast rational approximation to be later returned as info
RmFast = Rd;

% We return only one error: errLB if the user asked for a stopping
% criterion on the divided differences matrices; errvecMat if the
% research was exact

if exactLB
    errLB = errvecMat;
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

function b = bfun(k, zz, z, w)
% This function corresponds to the functions b_j(zz) =
% \frac{(zz-z_j)^{-1}}{\sum_{i=1}^d\frac{w_i}{zz-z_i}},
% where z and w are vectors of length d returned by the AAA algorithm and k
% is an index from 1 to d. Notice that b_k(z_k) = w_k^{-1}, so
% we force this interpolation.

% We encountered a really strange behaviour. If one computes b_k(z_0) (and
% only for z_0!), then he/she gets w_k^{-1}, instead of zero. This happens
% because MATLAB returns "Inf +NaNi" in the denominator, instead of
% "Inf+Infi". Consequently, num/den, i.e., 0/(Inf+NaNi) is NaN and not
% zero. Therefore the function thinks it has the form 0/0 and it
% interpolates. Why does this happen only for z_0? Not sure, but a strong
% suspicion lies in the fact that den = CC*w, where CC(1) is Inf and w(1)
% is real, while w(2:end) is complex (indeed, b_k(z_i) = 0, for i \neq 0 or
% k). Even weirder, we have CC*w = Inf +NaNi, but (when CC is a rwo vector)
% sum(CC.*w.') = Inf+Infi, as expected.
% We solved this issue with the command den(isnan(den)) = Inf.

zz = zz(:);
l = length(zz);
CC = 1./bsxfun(@minus,zz,z.'); % Cauchy matrix
den = CC*w; % denominator
% The next line solves a weird MATLAB behaviour, as explained at the
% beginning of the function
den(isnan(den)) = inf;
num = CC(:,k+1); % numerator, aka 1/(zz-z_k);
b = num./den;
% Now we force interpolation
ii = find(isnan(b));
for j=1:length(ii)
    b(ii(j)) = 1/w(find( zz(ii(j)) == z ));
end

end

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

function Rd = iRmatrixHandleLB(zz, RdStart, bLB, D)

Rd = RdStart(zz);
for k=1:length(D)
    Rd = Rd + bLB(k,zz)*D{k};
end

end

function Rd = iRmatrixHandleLBFast(zz, RdStart, bLB, D, DDense)

if issparse(D{1})
    bLBVec = zeros(length(D), 1);
    for k=1:length(D)
        bLBVec(k) = bLB(k,zz);
    end
    Rd = RdStart(zz) + DDense*bLBVec;
else
    Rd = RdStart(zz);
    for k=1:length(D)
        Rd = Rd + bLB(k,zz)*D{k};
    end
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
% Check whether Rd is NaN. If that is the case, it means zz = z_i for some
% i, therefore we force the interpolation.
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
    % Check whether Rd is NaN. If that is the case, it means zz = z_i for
    % some i, therefore we force the interpolation.
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



%% Parsing inputs
function Params = iParseInputs(F, Z, opts)

% Settings default parameters
dfVal = iDefaultValues(Z);
Params.phase2 = dfVal.phase2;
Params.tol1 = dfVal.tol1;
Params.tol2 = dfVal.tol2;
Params.Z2 = dfVal.Z2;
Params.dmax = dfVal.dmax;
Params.verbose = dfVal.verbose;
% End of default

if ~isa(F, 'struct')
    error('F must be a struct with fields "coeffs" and "fun".')
end

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
if isfield(opts, 'verbose')
    Params.verbose = opts.verbose;
end

% Check sparsity
Params.sparseFlag = 1;
for j = 1:length(F.coeffs)
    if ~issparse(F.coeffs{j})
        Params.sparseFlag = 0;
    end
end
%Convert Params.phase2 in the correspective easier parameters to parse
if isequal(Params.phase2,'exact')
    Params.exactSearch = 1;
    Params.LBref = 0;
    Params.exactLB  = 0;
elseif isequal(Params.phase2, 'relaxed')
    Params.exactSearch = 0.5;
    Params.LBref = 0;
    Params.exactLB  = 0;
elseif isequal(Params.phase2, 'LB')
    Params.exactSearch = 0;
    Params.LBref = 1;
    Params.exactLB  = 0;
elseif isequal(Params.phase2, 'exactLB')
    Params.exactSearch = 0;
    Params.LBref = 1;
    Params.exactLB  = 1;
elseif isequal(Params.phase2, 'relaxedLB')
    Params.exactSearch = 0;
    Params.LBref = 1;
    Params.exactLB  = 0.5;
elseif isequal(Params.phase2, '')
    Params.exactSearch = 0;
    Params.LBref = 0;
    Params.exactLB  = 0;
else
    error(['Input in opts.phase2 not recognized.'])
end

end

function DefaultValues = iDefaultValues(Z)

DefaultValues.tol1 = 1e-11;
DefaultValues.tol2 = 1e-13;
DefaultValues.Z2 = Z;
DefaultValues.phase2 = '';
DefaultValues.dmax = 20;
DefaultValues.verbose = 1;

end
