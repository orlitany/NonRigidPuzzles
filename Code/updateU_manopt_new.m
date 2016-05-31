function [u_hat,B] = updateU_manopt_new(previous_u,BMat,C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = numel(model.shape.X);

previous_u = vertcat(previous_u{:});
% previous_u = reshape(previous_u,[main_params.num_parts N]);

tmp = cellfun(@(X)X.model,corr_functions,'UniformOutput',false);
F = blkdiag(tmp{:});


CA = cellfun(@(X,Y)X*Y,C,A,'UniformOutput',false)';
CA = [CA{:}];

aggregate_mat = repmat(speye(N),[1 main_params.num_parts]);

area_mat = kron(eye(main_params.num_parts),model.shape.area_per_vertex');

regularizaion_mat = cellfun(@(X)(gradientNorm(model.shape,X)),rho,'UniformOutput', false);
regularizaion_mat = blkdiag(regularizaion_mat{:});

% cvx_begin quiet
% 
%     variable u(N*main_params.num_parts) nonnegative
%     minimize ( 0.5*square_pos (norm( CA - model.eigen_functions'*(aggregate_mat*diag(u))*F ,'fro'))  + 0.5*algo_params.lambda_s_model*u'*regularizaion_mat*u )
%     subject to
%         aggregate_mat*u == ones(N,1)    
%         area_mat*u >= 0.9*parts_areas;
% 
% cvx_end
% disp(['status for U update is: ' cvx_status])

manifold = euclideanfactory(main_params.num_parts*N);
problem.M = manifold;
% Phi = aggregate_mat'*model.eigen_functions;

problem.cost  =  @(x)fun_u_new(x,N,BMat,CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
problem.egrad = @(x)fun_du_new(x,N,BMat,CA,regularizaion_mat,area_mat,aggregate_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area,algo_params.lambda_ones)
                              

% Numerically check gradient consistency (optional).
% checkgradient(problem);
options.maxiter = 500;
u  = conjugategradient(problem,previous_u,options);
u_hat = num2cell(reshape(eta(u),[N main_params.num_parts]),1);

%- B
tmp = cellfun(@(X)size(X.model,2),corr_functions,'UniformOutput',false);
% B = mat2cell(model.eigen_functions'*aggregate_mat*diag(double(eta(u)))*F,main_params.num_eigen,[tmp{:}]);
B = mat2cell(reshape(BMat*eta(u),[size(CA,1) size(CA,2)]),main_params.num_eigen,[tmp{:}]);


end

