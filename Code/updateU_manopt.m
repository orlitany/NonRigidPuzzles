function [u_hat,B] = updateU_manopt(C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


N = numel(model.shape.X);


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

manifold = multinomialfactory(main_params.num_parts,N);
problem.M = manifold;
Phi = aggregate_mat'*model.eigen_functions;
problem.cost  =  @(x)fun_u(x,CA,Phi,F,regularizaion_mat,area_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area)
problem.egrad = @(x)fun_du(x,CA,Phi,F,regularizaion_mat,area_mat,parts_areas,algo_params.lambda_s_model,algo_params.lambda_area)
 
% Numerically check gradient consistency (optional).
% checkgradient(problem);
[u, xcost, info, options]  = conjugategradient(problem);

% Display some statistics.
figure;
for i=1:3,
    subplot(1,3,i);showshape(model.shape,u(i,:)');caxis([0 1]);colorbar
end
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');
u_hat = num2cell(reshape(u,[N main_params.num_parts]),1);

tmp = cellfun(@(X)size(X.model,2),corr_functions,'UniformOutput',false);

%- added a threshold of u
B = mat2cell(model.eigen_functions'*aggregate_mat*diag(double(u>0.5))*F,main_params.num_eigen,[tmp{:}]);



end

