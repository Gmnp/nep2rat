%This M-file tests different rational approximations of REPs and NEPs from NLEVP.

clear opts nepAcc nepDegree nepSteps

dmax = 60; %max degree
tol = 1e-13; %tolerance for phases 1 & 2
N = 10; %default problem size
nc = 300; %number of sample points for Sigma
nc2 = 50; %number of sample points on the countour of Sigma
tolnormest = tol/10;

% Z is a set of random points in the target set Sigma
% Z2 is a uniform discretisation of the contour of Sigma
% ZZ is the union of Z and Z2 minus (with no repetition)
useZZ = 1; %set to 1 to use ZZ for all the algs. (in particular phase 1)
useZ2 = 1; %set to 1 to use Z2 only in the refinement part (if not ZZ is used).

%Extract all reps and neps from NLEVP 4.1.
pep = nlevp('query','pep');
allProblems = nlevp('query','problems');
nep = setdiff(allProblems, pep);
nep = setdiff(nep, 'pillbox_cavity'); %Remove pillbox_cavity: too large
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nep = setdiff(nep, 'gun'); % Test gun on its own
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tol<1e-7 %remove the REPs - they are well approximated
   nep = setdiff(nep,'loaded_string');
   nep = setdiff(nep,'railtrack2_rep');
   nep = setdiff(nep,'railtrack_rep');
end
%nep= {'gun'} % When testing gun on its own
nb_test_pbs = length(nep);
fprintf('Number of test problems: %3.0f\n',nb_test_pbs)

% Option parameters for nep2rat.m
opts.dmax = dmax;
opts.tol1 = tol;
opts.tol2 = tol;
opts.tol  = tol;

%Arrays memory allocation
Z = zeros(nb_test_pbs,nc);
gam = zeros(nb_test_pbs,1); %centers
rad = zeros(nb_test_pbs,1); %radii
half_disc = gam;

% 1: set-valued AAA
% 2: weighted AAA
% 3: surrogate AAA
% 4: surrogate AAA with exact search
% 5: surrogate AAA with cyclic LB refinement
% 6: NLEIGS with cyclic LB procedure
alg_to_run = [2 3 4 5 6];

%Main loop
for kk = 1:nb_test_pbs
    opts.nF = [];
    switch nep{kk} %generate F

      case 'bent_beam' % temporary
            gam(kk) = 60;
            rad(kk) = 30;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'buckling_plate'  % The smallest poles are in k*pi/2, and in 4.70
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'canyon_particle'
            gam(kk) = -9e-2+1e-6i;
            rad(kk) = .1;
            half_disc(kk) = 1; %half disc domain
            stepSize = 1;
            [coeffs,fun,F] = nlevp(nep{kk}, stepSize);

      case 'clamped_beam_1d'
            gam(kk) = 0;
            rad(kk) = 10;
            [coeffs,fun,F] = nlevp(nep{kk}, 100);

      case 'distributed_delay1'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'fiber'
            gam(kk) = 0;
            rad(kk) = 2e-3;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'gun'
            gam(kk) = 62500;
            rad(kk) = 50000;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'hadeler'
            gam(kk) = -30;
            rad(kk) = 11.5;
            [coeffs,fun,F] = nlevp(nep{kk},200);

      case 'loaded_string'
            gam(kk) = 362;%14;
            rad(kk) = 358;%11;
            %KAPPA = 1; mass = 1; % pole is KAPPA/mass
            [coeffs,fun,F] = nlevp(nep{kk},100);%, N, KAPPA, mass);

      case 'nep1'
            gam(kk) = 0;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'nep2'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});

      case  'nep3'
            gam(kk) = 5i;
            rad(kk) = 2; % found 14 evals in this disc
            [coeffs,fun,F] = nlevp(nep{kk},N);

      case 'neuron_dde'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'photonic_crystal'
            % The poles of the default values are +-1.18-0.005i and +-1.26-0.01
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nep{kk}, N);

      case 'pillbox_small'
            gam(kk) = 0.08;
            rad(kk) = 0.05;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'railtrack2_rep'
            % railtrack2_rep has a pole in 0
            gam(kk) = 3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk}, N);

      case 'railtrack_rep'
            gam(kk) = -3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'sandwich_beam'
            gam(kk) = 7000;
            rad(kk) = 6900;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'schrodinger_abc'
            gam(kk) = -10;
            rad(kk) = 5;
            [coeffs,fun,F] = nlevp(nep{kk}, N);

      case 'square_root'
            gam(kk) = 10+50i;
            rad(kk) = 50;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'time_delay'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'time_delay2'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});

      case 'time_delay3'
            gam(kk) = 2;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nep{kk}, N, 5);

      otherwise
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk}, N);
    end

    Fsize(kk) = length(coeffs{1}); %record the size of each NEP
    %fprintf('Pb size %5d\n', Fsize(kk))
    fprintf('*******************************\n');
    fprintf('kk=%2d, Problem: %s, n =%4d\n',kk, nep{kk},Fsize(kk));
    fprintf('*******************************\n');

    %Generate set of points Z, Z2 and if needed ZZ = Z U Z2
    rng(0); %Fix the random number generator
    Z(kk,:) = rand(1,nc).*exp(rand(1,nc)*2*pi*1i);
    if half_disc(kk)
       negPoints = imag(Z(kk,:)) < 0;
       Z(kk,negPoints) = Z(kk, negPoints)';
       Z2 = gam(kk) + rad(kk)*exp(1i*linspace(0,pi,nc2)); % half circle
       Z2 = [Z2(2:end-1), linspace(-rad(kk), rad(kk), nc2)+gam(kk)];
    else
       Z2 = gam(kk) + rad(kk)*exp(1i*linspace(0,2*pi,2*nc2));
    end
    Z(kk,:) = Z(kk,:)*rad(kk)  + gam(kk); % shift to the correct points
    if useZZ  %merge Z and Z2
       ZZ = [Z(kk,:) Z2];
       %now look for repetitions and remove. We can't simply use "union"
       %because it is bugged
       ZRows = [real(ZZ)', imag(ZZ)'];
       Z1Rows = union(ZRows,ZRows,'rows');
       ZZ = Z1Rows(:,1)' + Z1Rows(:,2)'*1i;
    else
       ZZ = Z(kk,:);
    end
    if useZ2
       opts.Z2 = Z2;
    else
       opts.Z2 = ZZ;
    end

    %------------------------------------------
    %Now construct the different rational approx.
    %Use the Frobenius norm for (much) faster results
    %------------------------------------------
    alg = 0;
    
    %% Set valued AAA original
    disp('Set valued AAA original')
    alg = alg+1;
    algo_used{alg} = 'set valued AAA';
    [r, pol, res, zer, z, ff, w, errvec] = aaa_svOrig(fun, ZZ , 'tol', opts.tol2, 'mmax', opts.dmax+1);   
    Rm = @(z) iEvaluateRational(r, coeffs, z, issparse(coeffs{1}));
    nepDegree(kk,alg) = length(z)-1;
    nepSteps(kk,alg) = length(z)-1;
    %[nepAcc(kk,alg),normFZ] = computeApproxErr(F, Rm, ZZ, 'fro',Rm);
    [nepAcc(kk,alg),normFZ] = computeApproxErr(F, Rm, ZZ, 2, Rm, [], tolnormest);
    fprintf('||F||_S = %7.2e  (Sigma 2-norm)\n',normFZ);
  
    %% Weighted AAA
    disp('Weighted AAA')
    alg = alg+1;
    algo_used{alg} = 'weighted AAA';
    opts.phase2 = '';
    opts.phase1 = 'weighted';
    FWAAA.coeffs = coeffs;
    FWAAA.fun = fun;  
    [Am, Bm, Rm, info] = nep2rat(FWAAA, ZZ, opts); 
    nepDegree(kk,alg) = info.degree;
    nepSteps(kk,alg) = length(info.phase)-1;
    %[nepAcc(kk,alg),normFZ] = computeApproxErr(F, info.Rm, ZZ, 'fro', Rm);
    [nepAcc(kk,alg),normFZ]= computeApproxErr(F, Rm, ZZ, 2, Rm, [], tolnormest);
    fprintf('||F||_S = %7.2e  (Sigma 2-norm)\n',normFZ);

    %% Surrogate AAA
    disp('Surrogate AAA')
    alg = alg+1;
    algo_used{alg} ='surrogate';
    opts.phase1 = 'sur';
    opts.phase2 = '';
    [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);    
    nepDegree(kk,alg) = info.degree;
    nepSteps(kk,alg) = length(info.phase)-1;
    %nepAcc(kk,alg) = computeApproxErr(F, info.Rm, ZZ, 'fro', Rm);
    nepAcc(kk,alg) = computeApproxErr(F, Rm, ZZ, 2, Rm, normFZ,tolnormest);

    %% Surrogate AAA + Exact Search
    disp('Surrogate AAA+Exact')
    alg = alg+1;
    algo_used{alg} = 'surrogate + exact refinement';
    opts.phase2 = 'exact'; 
    [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
    nepDegree(kk,alg) = info.degree;
    nepSteps(kk,alg) = length(info.phase)-1;
    nepAcc(kk,alg) = computeApproxErr(F, Rm, ZZ, 2, Rm, normFZ,tolnormest);

    %% Surrogate AAA + LB
    disp('Surrogate AAA + LB')
    alg = alg+1;
    algo_used{alg} = 'surrogate+LB refinement';
    opts.phase2 = 'LB';
    opts.verbose = 0;
    opts.aaa = 0; 
    [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
    nepDegree(kk,alg) = info.degree;
    nepSteps(kk,alg) = length(info.phase)-1;
    %nepAcc(kk,alg) = computeApproxErr(F, info.Rm, ZZ, 'fro', Rm, normFZ);
    nepAcc(kk,alg) =  computeApproxErr(F, Rm, ZZ, 2, Rm, normFZ,tolnormest);

%%%%%%%%%%%
%NLEIGS can only be called after sur AAA+LB since it needs reordered poles.
% it also uses the lower bound on ||F||_Sigma from surrogateAAA (phase 1 only)
% stored in info.nF
    %% LB after we get ordered poles from surrogate,
    %% and cyclically repeat them
    disp('LB (NLEIGS) with poles from surrogate')
    alg = alg+1;
    algo_used{alg} = 'NLEIGS with surr AAA poles';
    opts.tol1 = tol;
    %ADD by GMNP
    if isempty(info.ordpol)
        %sandwich_beam returns a 0-th degree approximation if the tolerance
        %is large, so there are no poles and NLEIGS cannot run.
        setSpecialChar = 1;
    else
        setSpecialChar = 0;
        opts.Xi = info.ordpol(1:length(info.pol));
        opts.nF = info.nF;
    end
    if setSpecialChar
        nepDegree(kk,alg) = -1;
        nepAcc(kk,alg) = -1;
    else
        opts.cyclic = 1;
        opts.phase1 = 'LB';
        opts.phase2 = '';
        [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
        infoSL = info;
        nepDegree(kk,alg) = info.degree;
        warning off;
        %nepAcc(kk,alg) = computeApproxErr(F, Rm, ZZ, 'fro', Rm, normFZ);
        nepAcc(kk,alg) = computeApproxErr(F, Rm, ZZ, 2, Rm, normFZ,tolnormest);
    end
    nepSteps(kk,alg) = -1;

end

display('accuracy')
print_matrix(nepAcc,{'%7.0e','%7.0e','%7.1e','%7.1e','%7.1e','%7.1e'},[],7,1,1)
display('degree')
print_matrix(nepDegree,{'%3.0f','%3.0f','%3.0f','%3.0f','%3.0f','%3.0f'},[],7,1,1)
display('Steps')
print_matrix(nepSteps,{'%3.0f','%3.0f','%3.0f','%3.0f','%3.0f','%3.0f'},[],7,1,1)


