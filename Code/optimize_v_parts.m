function v = optimize_v_parts(M,N,corr_func,C,eta_u,mu1,mu2,opt,maxiter)

    % we abuse notation of this function to optimize for the part indicator
    % to do to we replace rolse:
    % A --> M.evecs'*M.S*G ; B --> C*N.evecs'*N.S ; 

    fprintf('Optimizing v...\n');

    A = M.evecs'*M.S*corr_func.model;
    B = C*N.evecs'*N.S;

    areas = full(diag(N.S));
    target_area = full( diag(M.S) )'*eta_u;
    x0 = N.evecs*C'*(M.evecs'*(M.S*eta_u));
%     x0 = M.evecs*((C*N.evecs')*(N.S*ones(N.n,1)));    
  
    manifold = euclideanfactory(N.n,1);

    problem = {};
    problem.M = manifold;

    vfunc = mumford_shah(N.VERT, N.TRIV, N.S);
    problem.cost =  @(v) vfunc.cost(A, B, corr_func.parts, v, ones(size(v)), target_area, areas, mu1, mu2, opt);
    problem.egrad = @(v) vfunc.grad(A, B, corr_func.parts, v, ones(size(v)), target_area, areas, mu1, mu2, opt,1);

%     figure
%     checkgradient(problem);

    options.maxiter = maxiter;%3e3;
    options.tolgradnorm = 1e-6;
    options.minstepsize = 1e-6;
    options.verbosity = 2;

    [v, cost, info, ~] = conjugategradient(problem, x0, options);

	fprintf('Optimize v, cost: %f\n', cost);
end
