function [Am, Bm, Rm, info] = nep2rat(F, Z, varargin)
%NEP2RAT Rational approximation to matrix-valued function and linearization.
% [Am, Bm, Rm, info] = NEP2RAT(F, Z) returns the rational
% approximation Rm of F and its linearization (Am, Bm).
% The exact algorithm depends on the form of F(z). User can gives optional
% parameters with the structure opts as [Am, Bm, Rm, info] = NEP2RAT(F,
% Z, opts)
% INPUT:
%  - F is the matrix valued function. It can either be a function_handle or
%  a structure that contains the matrix coefficients F.coeffs and the
%  scalar functions F.fun, as returned by the NLEVP library.
%  - Z is the set of points where the user wants the nonlinear function to
%  be approximated.
%  - Opts is an optional structure that contains all the optional
%  parameters for the code.
%     * .phase1: Tells the software which algorithm to use in the first
%     phase. If 'sur', it uses the surrogate AAA, if 'weighted', it uses
%     the weighted AAA, if 'LB', it uses the Leja-Bagby approximation. The
%     latter needs the parameter opts.Xi. The default is 'sur' for function
%     handle and 'weighted' for F being a struct.
%     * .phase2: Tells the software which algorithm to use in the second
%     phase. If 'exact' or 'relaxed', it uses an exact/relaxed search of
%     the best point for the algorithm used in the first phase. If 'LB', it
%     uses the Leja-Bagby refinement; if 'exactLB' or 'relaxedLB', it uses
%     an exact/relaxed search for Leja-Bagby; if '', it does not use any
%     refinement. The default value is 'LB' for F being a function_handle,
%     '' for F bein a struct.
%     * .dmax: Maximum degree of the approximation. Default is 20.
%     * .Z2: A subset of Z. If an exact or relaxed refinement is used, then
%     this refinement is done on  Z2. If Z2 is not a subset of Z, then the
%     software takes the union of Z and Z2 as the new Z.
%     * .tol1: The tolerance for the first phase. The default value is
%     1e-11.
%     * .tol2: The tolerance for the second phase. The default value is
%     1e-13.
%     * .verbose: If 1, it prints the warnings; if 2, makes the code more
%     verbose. Default is 0.
%     * .normalise: If 1, the algorithm works on the function G(z) =
%     F(z)/beta, where beta is an approximation of ||F||_Z = max_{z\in Z}
%     ||F(z)||_2. Default is 0.
%     * .Xi: If the first phase is 'LB', the algorithm needs a set of
%     possible poles Xi.
%     * .cyclic: If the first phase is 'LB', the poles can be either chosen
%     in the Leja-Bagby (0) way or extracted cyclically (1). Default is 0.
%     * .nF: If the first phase is 'LB', nF can be a lower bound of ||F||_Z
%     provided by the user. Default is 0.
%
% OUTPUT:
%  - Am and Bm are the pencil that linearize Bm.
%  - Rm is a matrix function handle that represents the rational
%  approximation.
%  - info returns the following information:
%     * .pol, .res, .zer: the poles, the residues and the zeros returned by
%     the AAA routine used.
%     * .z, .f, .w: the nodes, the function values and the weights
%     of the specific AAA routine used.
%     * .errvec: the errors of the phase 1 approximation.
%     * .errvec2: the errors of the phase 2 approximation.
%     * .phase2: a vector of the form [1,...,1,2,...,2], where the numbers
%     of 1s and 2s represents the steps in phase 1 and phase 2.
%     * .degree: The degree of the approximation.
%     * .msg: it contains any warning message.
%     * .ordPol, .zLB: if opts.phase2 is "LB*", then these contain the
%     all ordered poles and the sample points of the Leja-Bagby
%     approximation. The first parts of .z and .zLB coincide.

% References: S. Guettel, G. M. Negri Porzio, F. Tisseur, "Robust rational 
% approximation of nonlinear eigenvalue problems". Under review

%% Optional parameters
if nargin == 3
    opts = varargin{1};
else
    opts = struct();
end
Params = iParseInputs(F, Z, opts);

%% Setting the parameters
Z = Params.Z;
normalise = Params.normalise;
verbose = Params.verbose;
Xi = Params.Xi;
%% If Z and Xi have common elements, we cannot run the algorithm.
if ~isempty(intersect(Z,Xi))
    if verbose >= 1
        warning(['The algorithm does not work if Z and Xi have points ' ...
            'in common.']);
    end
    Am = [];
    Bm = [];
    info.msg = 'Z and Xi have points in common.';
    return
end

% End of optional parameters
if verbose >= 2
    fprintf("Approximating the nonlinear matrix-based function.\n")
end
%% Normalizing the function
if normalise
    [F, normFSigma] = iNormalize(F,Z);
    normFSigma;
end
%% Main core
if isa(F, 'function_handle')
    if isequal(Params.phase1,'sur')
        [Am, Bm, Rm, info] = hybrid_solver(F, Z, Params);
    elseif isequal(Params.phase1,'weighted')
        if verbose >= 1
            warning(['Algorithm cannot call weighted AAA if the matrix function',...
                ' is in function-handle form. Calling surrogate AAA instead.'])
        end
        [Am, Bm, Rm, info] = hybrid_solver(F, Z, Params);
    else % params.phase1 = LB
        [Am, Bm, Rm, info] = LejaBagbyRefinement(F, Z, Xi, Params);
    end
else % F is a struct
    Fstruct = F;
    F = @(z) iBuildFunctionHandle(F,z);
    if isequal(Params.phase1, 'sur')
        if verbose >=1
            warning('User called surrogate AAA even if the matrix function was in split form.')
        end
        [Am, Bm, Rm, info] = hybrid_solver(F, Z, Params);
    elseif isequal(Params.phase1, 'weighted')
        [Am, Bm, Rm, info] = weighted_aaaRef(Fstruct, Z, Params);
    else % Leja-Bagby refinement
        [Am, Bm, Rm, info] = LejaBagbyRefinement(F, Z, Xi, Params);
    end
end

%% Renormalize
if normalise
    Am = Am*normFSigma;
    Bm = Bm *normFSigma;
    Rm = @(z) Rm(z)*normFSigma;
    info.Rm = @(z) info.Rm(z)*normFSigma;
    info.r = @(z) info.r(z)*normFSigma;
    info.normFSigma = normFSigma;
end

end




%% Auxiliary functions

% If F is not a function_handle, iBuildFunctionHandle transforms it into
% one
function F = iBuildFunctionHandle(G, z)
coeffs = G.coeffs;
fun = G.fun;
funValues = fun(z);
d = length(coeffs);
F = funValues(1)*coeffs{1};
for j = 2:d
    F = F + funValues(j)*coeffs{j};
end

end

function [G, normFSigma] = iNormalize(F,Z)
% iNormalize(F,Z) computes an approximation normFSigma of ||F||_Z := max_{z\in Z}
% ||F(z)||_2 and return the matrix-valued function G(z) = F(z)/normFSigma.
% It disguishes the cases when F is given as a matrix-handle and a struct.
if isa(F, 'function_handle')
    n = size(F(Z(1)),1);
    state = rng(); rng('default');
    u = randn(n,1);
    u = u/norm(u);
    rng(state);
    normFSigma = 0;
    for jj =1:length(Z)
        normFSigma = max(normFSigma, norm(F(Z(jj))*u));
    end
    G = @(z) F(z)/normFSigma;
else % F is a struct
    sparseFlag = issparse(F.coeffs{1});
    n = size(F.coeffs{1},1);
    d = length(F.coeffs);
    state = rng();
    rng('default');
    u = randn(n,1);
    u = u/norm(u);
    rng(state);
    if sparseFlag
        uj = sparse(n,d);
    else
        uj = zeros(n,d);
    end
    for jj = 1:d
        uj(:,jj) = F.coeffs{jj}*u;
    end
    normFSigma = 0;
    for jj =1:length(Z)
        normFSigma = max(normFSigma, norm(uj*F.fun(Z(jj)).'));
    end
    % Normalized function
    G.coeffs = F.coeffs;
    G.fun = @(zz) F.fun(zz)/normFSigma;
end
end

%% Parsing inputs

% iParseInputs checks the auxiliary parameters
function Params = iParseInputs(F, Z, opts)

% Settings default parameters
dfVal = iDefaultValues(F, Z);
Params.phase1 = dfVal.phase1;
Params.phase2 = dfVal.phase2;
Params.tol1 = dfVal.tol1;
Params.tol2 = dfVal.tol2;
Params.Z2 = dfVal.Z2;
Params.Z = Z(:);
Params.dmax = dfVal.dmax;
Params.normalise = dfVal.normalise;
Params.Xi = dfVal.Xi;
Params.cyclic = dfVal.cyclic;
Params.verbose = dfVal.verbose;
% End of default


if isfield(opts, 'phase1')
    Params.phase1 = opts.phase1;
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
    Params.Z = union(Params.Z2, Params.Z, 'stable'); % We want Z being superset of Z2;
    % The following lines are needed because union does not always return a
    % vector with unique elements (MATLAB bug). Therefore we fix it
    % ourselves. As soon as we learn union is not bugged anymore, they can
    % be erased.
    equa = sparse(Params.Z-Params.Z.')==0;
    equa = equa-speye(size(equa)); % logical matrix where the elements are repeated;
    [ind, ~] = find(equa); % vector of the row indeces
    % due to the simmetry of equa, ind has a even length and contain both the
    %'original value' and the 'repeated one', so we erase only one copy;
    ind = ind(1:length(ind)/2);
    Params.Z(ind) = [];
    % end of lines to be cancelled
end
if isfield(opts, 'dmax')
    Params.dmax = opts.dmax;
end
if isfield(opts, 'normalise')
    Params.normalise = opts.normalise;
end
if isfield(opts, 'Xi') % This is necessary for LB in phase 1;
    Params.Xi = opts.Xi(:);
end
if isfield(opts, 'cyclic') % this is used only in the LB refinement, phase 1
    Params.cyclic = opts.cyclic;
end
if isfield(opts, 'verbose')
    Params.verbose = opts.verbose;
end
if isfield(opts, 'nF')
    Params.nF = opts.nF;
end
end

% Default parameters

function DefaultValues = iDefaultValues(F, Z)
if isa(F, 'function_handle')
    DefaultValues.phase1 = 'sur';
    DefaultValues.phase2 = 'LB';
elseif isa(F, 'struct')
    DefaultValues.phase1 = 'weighted';
    DefaultValues.phase2 = '';
else
    error(["Input not recognized. F must either be a struct with fields", ...
        " coeffs and fun, or a function_handle.\n"])
end

DefaultValues.Z2 = Z(:);
DefaultValues.tol1 = 1e-11;
DefaultValues.tol2 = 1e-13;
DefaultValues.dmax = 20;
DefaultValues.normalise = 0;
DefaultValues.nF = 0;
DefaultValues.Xi = [];
DefaultValues.cyclic = 0;
DefaultValues.verbose = 0;
end
