% This is an example where we approximate a scalar functions in different
% ways
clear opts
f = @(z) exp(z.^2);
% We test it for different values of alpha
for j = [-6:6]%3
    a = 10^j
    % We build the matrix valued function
    fun = @(z) [z*0+1, -z, f(z), -f(z)];
    realEvs= [0.5, 0.3];
    coeffs = {[realEvs(1) 0; 0 realEvs(2)], [1 0; 0 1], [a 1; 0 a], [a 10000; 0 a]};
    F.fun = fun;
    F.coeffs = coeffs;
    % random points in the target region, the unit disk
    nc = 500;
    Z = rand(1,nc).*exp(rand(1,nc)*2*pi*1i);
    % We want all the eigenvalues, so that we can show the spurious ones
    % are outside the target region
    opts.allEvs = 1;
    opts.verbose = 0;
    % The eigensolver itself
    [evs, evecs, resids, info] = mixedSolver(F, Z, opts);
    fprintf('Number of eigenvalues returned: %d\n', length(evs))
    fprintf('Number of eigenvalues returned inside the target region: %d\n', sum(abs(evs)<1))
    disp('Returned eigenvalues inside the target region: ')
    evsInRegion = evs(abs(evs)<1)   
end