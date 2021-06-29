function [evs, evecs, resids, info] = mixedSolver(F, Z, varargin)

% [Evs, Evecs, Resids, Info] = MIXED_SOLVER(F,Z,Opts) solves the nonlinear
% eigenvalue problem F(z)v = 0 through a rational approximation. The exact
% algorithm depends on the form of F(z).
% INPUT: 
%  - F is the matrix valued function. It can either be a function_handle or
%  a structure that contains the matrix coefficients F.coeffs and the
%  scalar functions F.fun, as returned by the NLEVP library.
%  - Z is the set of point where the user wants the nonlinear function to
%  be approximated.
%  - Opts is an optional structure that contains all the optional
%  parameters for the code. More on this below.
% OUTPUT:
%  - Evs, Evecs, Resids are respectively the eigenvalues, the eigenvectors
%  and the residuals of the problem.
%  - Info is a structure that contains additional information that mirrors
%  the fields of Opts.
%
% OPTS is a structure with the optional parameters:
%  - .evsThresh: The algorithm does not accept eigenvalues if their
%  residuals is larger than this value. Default is 1e-5.
%  - .nevs: If the input is sparse, nevs is the number of eigenvalues
%  returned by the algorithm. Default is 6.
%  - .aaa, .xi: If F is a function_handle, .aaa and .xi specify which
%  algorithm to rationally approximate F(z). See HELP HYBRID_SOLVER
%  for more info. Default is 0.5 and "EMPTY".
%  - .tol: The threshold for the rational AAA-approximation. Default
%  is 1e-13.
%  - .mmax: Max number of iteration of the AAA algorithm. Default is
%  20.
%  - .verbose: Fixing the verbosity of the code. Default value is 1
% [0,1].
 % - .Sigma: If the input is sparse, .Sigma is the nearest point where
 % eigslooks for eigenvalues. Default value is the centroid of all the
 % points in Z.
  


% Optional parameters
if nargin == 3
    opts = varargin{1};
else
    opts = struct();
end
Params = iParseInputs(Z, opts);
evsThresh = Params.evsThresh;
nevs = Params.nevs;
Sigma = Params.Sigma;
allEvs = Params.allEvs;
verbose = Params.verbose;

% End of optional parameters
if verbose
   fprintf("Approximating the nonlinear matrix-based function.\n") 
end

[Am, Bm, Rm, info] = nep2rat(F,Z,opts);

% This is the approximation of the norm ||F||_Z computed by nep2rat. We use
% it to compute the backward error of the eigenpairs

nF = info.nF;
% Retrieve eigenvalues and eigenvectors
% In the sparse case, if the number of eigenvalues was not given, we
% compute all the eigenvalues.

if isa(F, 'function_handle')
    n = size(F(Z(1)),1);
elseif  isa(F, 'struct')
     n = size(F.coeffs{1},1);
     Fstruct = F;
     F = @(z) iBuildFunctionHandle(F,z);
end

if strcmp(nevs, 'all')
      nevs = min(size(Am,1));
end

if verbose
   fprintf("Computing the eigenvalues.\n") 
end

%% Temporary %%
% Am = full(Am);
% Bm = full(Bm);
%% end       %%
if all(isfinite(Am), 'all') && all(isfinite(Bm), 'all')
    if issparse(Am)
        [evecs, evs] = eigs(Am, Bm, nevs, Sigma);
    else
        [evecs, evs] = eig(Am, Bm);
    end
    evs = diag(evs);
    evecs = evecs(1:n,:);
    % Normalize again the eigenvector in infty norm
    evecs = evecs./max(abs(evecs));
    % Check the eigenvalues with a large residuals 
    Nevs = length(evs);
    resids = zeros(Nevs, 1);
    badEvs = zeros(Nevs,1);
    for kk = 1:Nevs
        resids(kk) = norm(F(evs(kk))*evecs(:,kk))/(nF*norm(evecs(:,kk)));
        if isnan(resids(kk)) || (resids(kk) > evsThresh)
            badEvs(kk) = 1;
        end
    end
    % Remove the eigenvalues with a large residual if allEvs = 0
    if ~allEvs
        evs = evs(~badEvs);
        evecs = evecs(:, ~badEvs);
        resids = resids(~badEvs);
    else
        info.badEvs = badEvs;
        info.Rm2 = Rm;
        info.Am = Am;
        info.Bm = Bm;
    end
else
    evs = [];
    evecs = [];
    resids = [];
    info.Rm2 = Rm;
    info.Am = Am;
    info.Bm = Bm;
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

% iParseInputs checks the auxiliary parameters
function Params = iParseInputs(Z, opts)

% Settings default parameters
dfVal = iDefaultValues(Z);
Params.evsThresh = dfVal.evsThresh;
Params.nevs = dfVal.nevs;
Params.aaa = dfVal.aaa;
Params.xi = dfVal.xi;
Params.tol = dfVal.tol;
Params.dmax = dfVal.dmax;
Params.nevs = dfVal.nevs;
Params.verbose = dfVal.verbose;
Params.Sigma = dfVal.Sigma;
Params.allEvs = dfVal.allEvs;

% End of default


if isfield(opts, 'evsThresh')
    Params.evsThresh = opts.evsThresh;
end
if isfield(opts, 'aaa')
    Params.aaa = opts.aaa;
end
if isfield(opts, 'nevs')
Params.nevs = opts.nevs;
end
if isfield(opts, 'xi')
    Params.xi = opts.xi;
end
if isfield(opts, 'tol')
    Params.tol = opts.tol;
end
if isfield(opts, 'tol1')
    Params.tol1 = opts.tol1;
end
if isfield(opts, 'tol2')
    Params.tol2 = opts.tol2;
end
if isfield(opts, 'dmax')
    Params.dmax = opts.dmax;
end
if isfield(opts, 'nevs')
   Params.nevs = opts.nevs; 
end
if isfield(opts, 'verbose')
    Params.verbose = opts.verbose;
end
if isfield(opts, 'Sigma')
   Params.Sigma = opts.Sigma; 
end
if isfield(opts, 'allEvs')
   Params.allEvs = opts.allEvs; 
end


end

% Default parameters

function DefaultValues = iDefaultValues(Z)

DefaultValues.evsThresh = 1e-5;
DefaultValues.nevs = 6;
DefaultValues.aaa = 0.5;
DefaultValues.xi = 0.5;
DefaultValues.tol = 1e-13;
DefaultValues.tol1 = 1e-11;
DefaultValues.tol2 = 1e-13;
DefaultValues.dmax = 20;
DefaultValues.nevs = 'all';
DefaultValues.verbose = 1;
DefaultValues.Sigma = sum(Z)/length(Z);
DefaultValues.allEvs = 0;

end
