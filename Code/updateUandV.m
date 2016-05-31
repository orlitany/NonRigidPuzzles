function [x] = updateUandV(previous_x,SSP_blk,BMat,C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = numel(model.shape.X);
K = size(previous_x{1},1);

previous_alpha = vertcat(previous_x{:});

for prt=1:main_params.num_parts
    Q{prt} = precomputeB_fast(parts{prt}.eigen_functions'*parts{prt}.M,corr_functions{prt}.parts);
end
CAMat = blkdiag(Q{:});

Q1 = [-BMat CAMat];


aggregate_mat = repmat(speye(N),[1 main_params.num_parts]);

area_mat = kron(eye(main_params.num_parts),model.shape.area_per_vertex');

regularizaion_mat = cellfun(@(X)(gradientNorm(model.shape,X)),rho,'UniformOutput', false);
regularizaion_mat = blkdiag(regularizaion_mat{:});

manifold = euclideanfactory(main_params.num_parts*K);
problem.M = manifold;
% Phi = aggregate_mat'*model.eigen_functions;
BMat_t_CA =   BMat'*CA(:);
problem.cost  =  @(x)fun_u_SSP(x,N,SSP_blk,BMat,CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
problem.egrad = @(x)fun_du_SSP_fast(x,N,SSP_blk,BMat,BMat_t_CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
                              

% Numerically check gradient consistency (optional).
% checkgradient(problem);
options.maxiter = algo_params.num_iter.manopt_maxiter;
alpha  = conjugategradient(problem,previous_alpha,options);
alpha_hat = num2cell(reshape(alpha,[K main_params.num_parts]),1);

end

