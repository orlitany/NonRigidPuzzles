function [alpha_hat] = updateU_manopt_SSP(previous_alpha,SSP_blk,BMat,C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = numel(model.shape.X);
K = main_params.num_subspace;

if ~isempty(previous_alpha), 
    previous_alpha = vertcat(previous_alpha{:});
end
% previous_u = reshape(previous_u,[main_params.num_parts N]);

% tmp = cellfun(@(X)X.model,corr_functions,'UniformOutput',false);
% F = blkdiag(tmp{:});


CA = cellfun(@(X,Y)X*Y,C,A,'UniformOutput',false)';
CA = [CA{:}];

aggregate_mat = repmat(speye(N),[1 main_params.num_parts]);

area_mat = kron(eye(main_params.num_parts - main_params.is_missing_part),model.shape.area_per_vertex');
if main_params.is_missing_part,
    area_mat = [area_mat zeros(size(area_mat,1),N)];
end

regularizaion_mat = cellfun(@(X)(gradientNorm(model.shape,X)),rho,'UniformOutput', false);
regularizaion_mat = blkdiag(regularizaion_mat{:});

manifold = euclideanfactory(main_params.num_parts*K);
problem.M = manifold;
% Phi = aggregate_mat'*model.eigen_functions;
BMat_t_CA =   BMat'*CA(:);
SSP_blk = sparse(SSP_blk);
BMat = sparse(BMat);
problem.cost  =  @(x)fun_u_SSP(x,N,SSP_blk,BMat,CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
problem.egrad = @(x)fun_du_SSP_fast(x,N,SSP_blk,BMat,BMat_t_CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
                              

% Numerically check gradient consistency (optional).
% checkgradient(problem);
options.maxiter = algo_params.num_iter.manopt_maxiter;
% profile on 
alpha  = conjugategradient(problem,previous_alpha,options);
% profile viewer
alpha_hat = num2cell(reshape(alpha,[K main_params.num_parts]),1);

end

