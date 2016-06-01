function u = optimize_u_multi_parts(u,M,N,corr_func,C,target_area,mu1,mu2,lambda_ones,opt,maxiter)

    num_parts = numel(u);
    u = cell2mat(u);
    G = cellfun(@(X)X.model,corr_func,'UniformOutput',false);
    fprintf('Optimizing v...\n');
    
    A = cellfun(@(CC,NN,F)CC*NN.evecs'*NN.S*F.parts,C,N,corr_func,'UniformOutput',false);
    B = M.evecs'*M.S;

    areas = full(diag(M.S));
%     target_area = full( sum(diag(N.S)) );

%     x0 = M.evecs*((C*N.evecs')*(N.S*ones(N.n,1)));
    %x0 = ones(M.n,1);
    
    manifold = euclideanfactory(M.n,num_parts);
    

    problem = {};
    problem.M = manifold;

    vfunc = mumford_shah(M.VERT, M.TRIV, M.S);
%     problem.cost =  @(v) vfunc.cost(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);
%     problem.egrad = @(v) vfunc.grad(A, B, G, v, ones(size(v)), target_area, areas, mu1, mu2, opt);
    problem.cost = @(x)multiparts_u(@vfunc.cost, A, B, G, x, target_area, areas, mu1, mu2,lambda_ones, opt);
    problem.egrad = @(x)multiparts_du(@vfunc.grad, A, B, G, x, target_area, areas, mu1, mu2, lambda_ones,opt);

%     figure
%     checkgradient(problem);

    options.maxiter = maxiter;%5e2;%3e3;
    options.tolgradnorm = 1e-6;
    options.minstepsize = 1e-6;
    options.verbosity = 2;

    [u, cost, info, ~] = conjugategradient(problem, u, options);
    u = num2cell(u,1);
	fprintf('Optimize v, cost: %f\n', cost);
end
