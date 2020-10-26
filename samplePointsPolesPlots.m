nepPlot = {'bent_beam', 'buckling_plate', 'gun', 'nep1', 'schrodinger_abc'};
nepPlot = setdiff(nepPlot, 'bent_beam');
dmax = 60; %max degree
tol = 1e-10; %tolerance for phases 1 & 2
N = 10; %default problem size
%nc = 1000; %number of sample points for Sigma
nc = 300; %number of sample points for Sigma
nc2 = 100;
useZZ = 1;
useZ2 = 1;
nb_test_pbs = length(nepPlot);
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

% Options for the plot
beta = 3.5; % stretch for the axes. We can probably set a specific one for each problem
nRows = ceil(sqrt(nb_test_pbs));
nCols = ceil(nb_test_pbs/nRows);
mkSizeSigma = 5;
mkPts = 7;
hold off
for kk = 1:nb_test_pbs
    switch nepPlot{kk} %generate F

      case 'bent_beam' % temporary
            gam(kk) = 60;
            rad(kk) = 30;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'buckling_plate'  % The smallest poles are in k*pi/2, and in 4.70
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'canyon_particle'
            gam(kk) = -9e-2+1e-6i;
            rad(kk) = .1;
            half_disc(kk) = 1; %half disc domain
            stepSize = 1;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, stepSize);

      case 'clamped_beam_1d'
            gam(kk) = 0;
            rad(kk) = 10;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, 100);

      case 'distributed_delay1'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'fiber'
            gam(kk) = 0;
            rad(kk) = 2e-3;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'gun'
            gam(kk) = 62500;
            rad(kk) = 50000;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'hadeler'
            gam(kk) = -30;
            rad(kk) = 11.5;
            [coeffs,fun,F] = nlevp(nepPlot{kk},200);

      case 'loaded_string'
            gam(kk) = 362;%14;
            rad(kk) = 358;%11;
            %KAPPA = 1; mass = 1; % pole is KAPPA/mass
            [coeffs,fun,F] = nlevp(nepPlot{kk},100);%, N, KAPPA, mass);

      case 'nep1'
            gam(kk) = 0;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'nep2'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case  'nep3'
            gam(kk) = 5i;
            rad(kk) = 2; % found 14 evals in this disc
            [coeffs,fun,F] = nlevp(nepPlot{kk},N);

      case 'neuron_dde'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'photonic_crystal'
            % The poles of the default values are +-1.18-0.005i and +-1.26-0.01
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, N);

      case 'pillbox_small'
            gam(kk) = 0.08;
            rad(kk) = 0.05;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'railtrack2_rep'
            % railtrack2_rep has a pole in 0, so we shifted it in 5 (?)
            gam(kk) = 3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, N);

      case 'railtrack_rep'
            gam(kk) = -3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'sandwich_beam'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'schrodinger_abc'
            gam(kk) = -10;
            rad(kk) = 5;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, N);

      case 'square_root'
            gam(kk) = 10+50i;
            rad(kk) = 50;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'time_delay'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'time_delay2'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nepPlot{kk});

      case 'time_delay3'
            gam(kk) = 2;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, N, 5);

      otherwise
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nepPlot{kk}, N);
    end
   
    Fsize(kk) = length(coeffs{1}); %record the size of each NEP
    %fprintf('Pb size %5d\n', Fsize(kk))
    fprintf('*******************************\n');
    fprintf('kk=%2d, Problem: %s, n =%4d\n',kk, nepPlot{kk},Fsize(kk));
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
       %now look for repetitions and remove
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

      alg = 0;

    %% Surrogate AAA
    disp('Surrogate AAA')
    alg = alg+1;
    algo_used{alg} ='surrogate';
    opts.phase1 = 'sur';
    opts.phase2 = '';
    [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
    subplot(nRows,nCols, kk)
    hold off
    plot(ZZ, '.b', 'MarkerSize',mkSizeSigma);
    hold on
    plot(info.z, 'or', 'MarkerSize',mkPts)
    plot(info.pol, 'xr', 'MarkerSize',mkPts)
    
    
%     %% Set valued AAA original
%     disp('Set valued AAA original')
%     alg = alg+1;
%     algo_used{alg} = 'set valued AAA';
%     [r, pol, res, zer, z, ff, w, errvec] = aaa_svOrig(fun, ZZ , 'tol', opts.tol2, 'mmax', opts.dmax+1);
%     plot(z, 'o');
%     plot(pol, 'x');
    
    %% Weighted AAA
    disp('Weighted AAA')
    alg = alg+1;
    algo_used{alg} = 'weighted AAA';
    opts.phase2 = '';
    opts.phase1 = 'weighted';
    FWAAA.coeffs = coeffs;
    FWAAA.fun = fun;
    [Am, Bm, Rm, info] = nep2rat(FWAAA, ZZ, opts);
    plot(info.z, 'sg', 'MarkerSize',mkPts);
    plot(info.pol, '+g', 'MarkerSize',mkPts);
    
    
    %% Surrogate AAA + LB
    disp('Surrogate AAA + LB')
    alg = alg+1;
    algo_used{alg} = 'surrogate+LB refinement';
    opts.phase2 = 'LB';
    [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
    plot(info.zLB(length(info.z)+1:end), 'dk', 'MarkerSize',mkPts);
    
%     
%     %% LB after we get ordered poles from surrogate,
%     %% and cyclically repeat them
%     disp('LB (NLEIGS) with poles from surrogate')
%     alg = alg+1;
%     algo_used{alg} = 'NLEIGS with surr AAA poles';
%     %ADD by GMNP
%     if ~isempty(info.msg) || isempty(info.ordpol)
%         %sandwich_beam returns a 0-th degree approximation if the tolerance
%         %is large, so there are no poles and NLEIGS cannot run.
%         setSpecialChar = 1;
%     else
%         setSpecialChar = 0;
%         opts.Xi = info.ordpol(1:length(info.pol));
%     end
%     if ~setSpecialChar     
%         opts.cyclic = 1;
%         opts.phase1 = 'LB';
%         opts.phase2 = '';
%         [Am, Bm, Rm, info] = nep2rat(F, ZZ, opts);
%         plot(info.zLB, 'o')
%     end
    % nice title in tex
    mytitle = replace(nepPlot{kk}, '_', '\_');
    title(mytitle, 'Interpreter', 'latex')
    switch nepPlot{kk}
        
        case 'bent_beam'
               axis([-80 95 -12 35])
    
        case 'buckling_plate'
                axis([-15 35 -10 10])


        case 'gun'
               axis([-90000  117000 -50000 55000])
        case 'nep1'
                 axis([-5 5 -5 5])

        case 'schrodinger_abc'
                axis([-17 10 -25 25])
        otherwise
            axis([real(gam(kk))-rad(kk)*beta real(gam(kk))+rad(kk)*beta imag(gam(kk))-rad(kk)*beta imag(gam(kk))+rad(kk)*beta])

    end
end
%      legend({'Set $\Sigma$', '$\sigma_i$ for Surrogate AAA', ' $\xi_i$ for Surrogate AAA',  '$\sigma_i$ for set-valued AAA', ...
%     '$\xi_i$ for set-valued AAA',  '$\sigma_i$ for weighted AAA', '$\xi_i$ for weighted AAA', '$\sigma_i$ added by LB refinement', '$\sigma_i$ for SLEPc'}, 'Interpreter', 'latex')
  figure(1)  
 legend({'Set $\Sigma$', '$\sigma_i$ for surr AAA', '$\xi_i$ for surr AAA', ...
   '$\sigma_i$ for wght AAA', '$\xi_i$ for wght AAA', '$\sigma_i$ added by LB'}, 'Interpreter', 'latex', 'Location','SouthWest')
    
