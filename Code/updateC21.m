function C = updateC21(C,A,B,W,d,D,main_params,algo_params)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
mu_1 = algo_params.mu_1;
mu_2 = algo_params.mu_2;

k = size(C{1},1);
options.maxiter = algo_params.num_iter.manopt_maxiter_C;
options.tolgradnorm = 1e-06;
options.minstepsize = 1e-06;
parfor prt = 1:main_params.num_parts    
    problem = {};
    problem.M = euclideanfactory(k,k);
    problem.cost  =  @(C)(sum(sqrt(sum((C*A{prt} - B{prt}).^2,1))) + mu_1*norm(C.*W{prt},'fro')^2 ...
        +mu_2 * (norm(C'*C,'fro')^2 - sum(diag(C'*C).^2) + sum((diag(C'*C) - d{prt}').^2) ));
    problem.egrad = @(C)(norm_21_gradient(C,A{prt},B{prt}) + mu_1*2*W{prt}.^2.*C + mu_2 * 4*(C*C'*C - C.*D{prt} ));
    figure;checkgradient(problem);
    C{prt} = conjugategradient(problem,C{prt},options);    
end
