function [alpha_hat] = updateU21(previous_alpha,SSP,C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(previous_alpha)
    previous_alpha = cell2mat(previous_alpha);
end

N = numel(model.shape.X);
K = main_params.num_subspace;

CA = cellfun(@(X,Y)X*Y,C,A,'UniformOutput', false);

% aggregate_mat = repmat(speye(N),[1 main_params.num_parts]);
% 
% area_mat = kron(eye(main_params.num_parts - main_params.is_missing_part),model.shape.area_per_vertex');
% if main_params.is_missing_part,
%     area_mat = [area_mat zeros(size(area_mat,1),N)];
% end
%
area_mat = model.shape.area_per_vertex';
regularizaion_mat = cellfun(@(X)(gradientNorm(model.shape,X)),rho,'UniformOutput', false);
% regularizaion_mat = blkdiag(regularizaion_mat{:});

problem.M = euclideanfactory(K,main_params.num_parts);                   
problem.cost  =  @(x)fun_u21(x,SSP,model.eigen_functions'*model.M,CA,corr_functions{1}.model,regularizaion_mat,algo_params.lambda_s_model,area_mat,parts_areas,algo_params.lambda_area);
problem.egrad = @(x)fun_du21(x,SSP,model.eigen_functions'*model.M,CA,corr_functions{1}.model,regularizaion_mat,algo_params.lambda_s_model,area_mat,parts_areas,algo_params.lambda_area);
                              

%- Numerically check gradient consistency (optional).
% checkgradient(problem)
options.maxiter = algo_params.num_iter.manopt_maxiter;
% profile on 
alpha  = conjugategradient(problem,previous_alpha,options);
% profile viewer
alpha_hat = num2cell(reshape(alpha,[K main_params.num_parts]),1);

end

