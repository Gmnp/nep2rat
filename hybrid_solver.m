function [Am, Bm, Rm, info] = hybrid_solver(F, Z, varargin)
% HYBRID_SOLVER    Rational approximant for nonlinear eigenvalue problems
% [Am,Bm,Rm,info] = HYBRID_SOLVER(F,Z,OPTS) is a rational approximant for
% nonlinear matrix-valued function that exploits both the aaa method and
% the Leja-Bagby approach. This is an AUXILIARY function called by NEP2RAT,
% it should not be used on its own. See HELP NEP2RAT for more information.

% INPUTS: F: matrix-valued function handle to approximate.
%         Z: l-by-1 vector of points we use to approximate F.
%         OPTS: struct that contains the optional value that tells
%         the algorithm how to approximate F. More on that below.
%
%
%         - opts.tol1: First tolerance for the surrogate aaa algorithm.
%           Default is 1e-11;
%         - opts.tol2: Tolerance for Leja-Bagby refinement. Default value
%           is 1e-13.
%         - opts.exactSearch: If opts.aaa 0.5, then surrogate AAA
%           terminates by doing an exactSearch of the best matrix point
%           when it is equal to 1. If exactSearch is 0.5, then it does a
%           relaxed search. If 0, it performs the normal Leja-Bagby
%           refinement. Default is 0.
%         - opts.verbose: if 1, then the software prints more information.
%           Default is 0.


if nargin == 3
    opts = varargin{1};
else
    opts = struct();
end
Params = iParseInputs(F, Z, opts);
sparseFlag = issparse(F(Z(1)));

%Setting optional parameters
tol1 = Params.tol1;
tol2 = Params.tol2;
Z2 = Params.Z2;
dmax = Params.dmax;
verbose = Params.verbose;
exactSearch = Params.exactSearch;
exactLB = Params.exactLB;
LBref = Params.LBref;
%End of optional parameters

%Cleaning the input
Z = Z(:);
l = length(Z);
l2 = length(Z2);
[m,n] = size(F(1));
ff = zeros(l,1);

if verbose >= 2
    fprintf('Using hybrid_solver.\n')
end

state = rng(); rng('default');
u = randn(m,1);
u = u/norm(u);
v = randn(n,1);
v = v/norm(v);
rng(state);

FZ2 = cell(1, l2);
% We save FZ2 = F(Z2(i)) if we are using it in the exact Search or in the
% exact LB
if exactSearch || exactLB
    FZ2 = arrayfun(F,Z2,'UniformOutput',false);
end

nF = 0;
for j = 1:l
    temp = F(Z(j))*v;
    ff(j) = u'*temp;
    nF = max(nF,norm(temp));
end
%nF is a lower bound on the sigma-2-norm of F. 

% Calling AAA on the surrogate function. Notice that the input is dmax+1,
% because aaaSearch.m wants the maximum number of points.
[~,pol,res,zer,z,f,w,errvec, errvecMat] = aaaSearch(F,ff,Z,FZ2,Z2,...
    tol1,tol2*nF,dmax+1,exactSearch,verbose);


% Setting the returned information
if isempty(pol)
    pol = inf;
end
info.pol = pol;
info.res = res;
info.zer = zer;
info.z = z;
info.f = f;
info.w = w;
info.errvec = errvec;
info.errvec2 = errvecMat/nF;
phase = [ones(length(errvec),1); 2*ones(length(errvecMat),1)];
info.phase = phase;
info.Z2 = Z2;
% Degree of the approximation
d = length(z) - 1;
info.degree = d;
info.msg = '';
info.nF = nF;

% Usually dPol = d, but in this we take care of F being a polynomial
dPol = length(pol);

% We want to return a warning if we reached the max degree of
% approximation. Due to the cleaning up of the Froissart doublets, the
% degree could be lower than the max degree, even if we reached dmax.
% Therefore the criterion with the vector phase.
if length(info.phase) == dmax+1
    if isempty(find(phase == 2, 1))
        if verbose >= 1
            warning(['The algorithm reached the maximum number of iterations', ...
                ' during the surrogate AAA portion'])
        end
        info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the surrogate AAA portion'];
    else
        if verbose >= 1
            warning(['The algorithm reached the maximum number of iterations', ...
                ' during the surrogate AAA with exact/relaxed search portion'])
        end
        info.msg = ['The algorithm reached the maximum number of iterations', ...
            ' during the surrogate AAA with exact/relaxed search portion'];
    end
end



% matrices F(z_i), to build the rational approximation in barycentric form.
% They are store in a cell of d+1 m-by-n matrices
Fz = arrayfun(F,z,'UniformOutput',false);

% This is a constant to be used in the second loop stopping criterion. If
% we are doing an exact search in the Leja-Bagby refinement as the second
% phase, then we compute the norms of all the points in the Z2 vector.
if exactLB
    maxF = max(cellfun(@(z) norm(z, 'fro'), FZ2));
else
    maxF = max(cellfun(@(z) norm(z, 'fro'), Fz));
end
maxF = max(maxF/sqrt(n),nF); % Sharper lower bound on Sigma-2-norm of F

% Build the matrix-handle of the approximation in the fast way: if the
% matrix is sparse, then it is densified. It is used for testing
if sparseFlag
    indDense = find(Fz{1});
    FzDense = zeros(length(indDense), length(Fz));
    for j = 1:length(Fz)
        FzDense(:,j) = Fz{j}(indDense);
    end
else
    FzDense = []; % tokens for iRmatrixHandleFast
    indDense = []; % idem
end
Rd = @(lambda) iRmatrixHandleFast(lambda,z,Fz,w,FzDense,indDense);
info.Rm = Rd;

% If we have done an exact search or we don't do the Leja-Bagby refinement,
% or we already reached dmax, it means we have finished, so we need to
% build Am, Bm and return
if (exactSearch || ~LBref || d == dmax)
    Rm = @(lambda) iRmatrixHandle(lambda,z,Fz, w);
    % iLinearize takes care of the linearization, both in the simple
    % surrogate_aaa case and when the linearization is mixed. It
    % distingueshes between the two using d and dmax: when they are equal,
    % it means it did not enter the Leja-Bagby refinement, this is why we
    % are giving two "ds" in input. In addition, it would need the vector
    % of beta and of poles xi with the mixed linearization, so here they
    % are empty.
    if verbose >= 2
        fprintf('Linearizing the rational approximation.\n')
    end
    
    [Am, Bm] = iLinearize(z,w,Fz,[],[],m,n,d,d, sparseFlag);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we start Leja-Bagby refinement %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose >=2
    fprintf('Starting the Leja-Bagby refinement.\n')
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

[Rm, z, w, Fz, beta, xi, errLB, RmFast] = iLejaBagbyRefinement(F, FZ2, ...
    Rd, Z, Z2, z, w, pol, Fz, dmax, maxF, tol2, sparseFlag, exactLB);

info.Rm = RmFast;

%Linearization
if verbose >=2
    fprintf('Linearizing the rational approximation.\n')
end
[Am, Bm] = iLinearize(z,w,Fz,beta,xi,m,n,d,dmax, sparseFlag);
info.ordpol = xi;
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


function [Rm, z, w, Fz, beta, xi, errLB, RmFast] = iLejaBagbyRefinement...
    (F, FZ2, Rd, Z, Z2, z, w, pol, Fz, dmax, maxF, tol2, mysparse, exactLB)

l = length(Z);
l2 = length(Z2);
d = length(z)-1;
dPol = length(pol);

bPol = ones(dPol,1);
bZ = ones(l,1);
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
            %them out.
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
    n = length(Fsigma); 
    D{j-d} = (Fsigma - RdSparse(z(j+1)) )/bZM(ind);
    % Tighter stopping criterion.
    nFF = norm(Fsigma, 'fro')/sqrt(n);
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
      % errLB(j-d) = [];
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

function [Am, Bm] = iLinearize(z, w, Fz , beta, xi, m,n, d,dmax, mysparse)
% This auxiliary function linearize the rational approximation. We need two
% distinctions. The first one is caused by sparsity, because in the full
% case, Fz is a 3D array containing the matrices  of the evaluations F(z_i)
% and the divided differences matrices D_i, while in the sparse case, Fz is
% a cell. The second one is on the degree. If d=dmax, it means that we did
% not enter the Leja-Bagby refinement, therefore the linearization is a bit
% simpler.

% Another way we did not enter the Leja-Bagby refinement is if
% beta is empty. In order to avoid worst corrections in other part of the
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
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + firstBlkRow(:,end-n+1:end);
        % Build the diagonal of Am
        aux = diag(z(1:end-1)) + diag(z(1:end-2), -1);
        Am = kron(aux,I);
        % Build the first row of Am. We already have F(z_i)*w_i
        Am(1:m,:) = cell2mat(FFw(1:end-1));
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + z(end-1)*FFw{end}/z(end);
    else
        % we entered the refinement cycle, so the linearization is different
        beta = beta(:);
        N = length(z)-1; % k+l in the notes.pdf
        I = speye(m,n);
        % Build the diagonal and subdiagonal part of Bm
        aux = diag([ones(d, 1); 1./beta(1:end-1)], -1);
        aux = aux -diag(ones(N,1));
        Bm = kron(aux,I);
        % Build the first row of Bm
        FFw = Fz.';
        
        for j = 1:d+1
            FFw{j} = Fz{j}*w(j); %first d+1 term: %F(z_i)*w_i
        end
        firstBlkRow = cell2mat(FFw)/xi(end);
        Bm(1:m,:) = firstBlkRow(:,1:end-n);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + firstBlkRow(:,end-n+1:end)/beta(end);
        % Build the diagonal of Am
        aux1 = z(d+1:end-2)./beta(1:end-1);
        aux = diag([z(1:d); aux1], -1);
        aux = aux - diag([z(1:d+1); xi(d+1:end-1)]);
        Am = kron(aux,I);
        % Build the first row of Am
        Am(1:m,:) = cell2mat(FFw(1:end-1));
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + Fz{end}*z(end-1)/(beta(end)*xi(end));
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
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + z(end-1)*FFw(:,:,end)/z(end);
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
        FFw(:,:,1:d+1) = bsxfun(@times, Fz(:,:,1:d+1), reshape(w, [1,1,d+1])); %F(z_i)*w_i
        firstBlkRow = bsxfun(@rdivide, FFw, xi(end));
        Bm(1:m,:) = reshape(firstBlkRow(:,:,1:end-1), [m,n*N]);
        Bm(1:m, end-n+1:end) = Bm(1:m, end-n+1:end) + firstBlkRow(:,:,end)/beta(end);
        
        % Build the diagonal of Am
        aux1 = z(d+1:end-2)./beta(1:end-1);
        aux = diag([z(1:d); aux1], -1);
        aux = aux - diag([z(1:d+1); xi(d+1:end-1)]);
        Am = kron(aux,I);
        % Build the first row of Am
        FFw = Fz;
        FFw(:,:,1:d+1) = bsxfun(@times, Fz(:,:,1:d+1), reshape(w, [1,1,d+1])); %F(z_i)*w_i
        Am(1:m,:) = reshape(FFw(:,:,1:end-1), [m,n*N]);
        Am(1:m, end-n+1:end) = Am(1:m, end-n+1:end) + Fz(:,:,end)*z(end-1)/(beta(end)*xi(end));
    end
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
%CC = 1./(zz-z.');
den = CC*w; % denominator
%The next line corrects a weird behaviour of MATLAB, as explained at the
%beginning of the function
den(isnan(den)) = inf;
num = CC(:,k+1); % numerator, aka 1/(zz-z_k);
b = num./den;
% Now we force interpolation
ii = find(isnan(b));
for j=1:length(ii)
    b(ii(j)) = 1/w(find( zz(ii(j)) == z ));
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

%% Parse inputs
% iParseInputs checks the auxiliary parameters
function Params = iParseInputs(F, Z, opts)

% Settings default parameters
dfVal = iDefaultValues(F, Z);
Params.phase2 = dfVal.phase2;
Params.tol1 = dfVal.tol1;
Params.tol2 = dfVal.tol2;
Params.Z2 = dfVal.Z2;
Params.dmax = dfVal.dmax;
Params.verbose = dfVal.verbose;
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
if isfield(opts, 'verbose')
    Params.verbose = opts.verbose;
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

% Default parameters

function DefaultValues = iDefaultValues(F, Z)
if ~isa(F, 'function_handle')
    error(["Input not recognized. F must either be a function_handle.\n"])
end
DefaultValues.phase2 = 'LB';
DefaultValues.tol1 = 1e-11;
DefaultValues.tol2 = 1e-13;
DefaultValues.Z2 = Z;
DefaultValues.dmax = 20;
DefaultValues.verbose = 1;
end
