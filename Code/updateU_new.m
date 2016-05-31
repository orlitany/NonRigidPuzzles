function [u_hat,B] = updateU_new(BMat,C,A,rho,parts_areas,model,corr_functions,main_params,algo_params)
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

% D=gradientNorm(model.shape,[])
% save('toMichael.mat','D','aggregate_mat','CA','F','Phi','area_mat','parts_areas')

cvx_begin quiet

    variable u(N*main_params.num_parts) nonnegative
%     minimize ( sum( norms( CA - reshape(BMat*u,[size(CA,1) size(CA,2)]) ,2,1)  ) + 0.5*algo_params.lambda_s_model*u'*regularizaion_mat*u )
    minimize ( 0.5*square_pos (norm( CA(:) - BMat*u ,'fro'))  + 0.5*algo_params.lambda_s_model*u'*regularizaion_mat*u  + 0.5*algo_params.lambda_area*square_pos(norm(area_mat*u - parts_areas,'fro')) );
    subject to
        aggregate_mat*u == ones(N,1)    
        

cvx_end
disp(['status for U update is: ' cvx_status])
u_hat = num2cell(reshape(u,[N main_params.num_parts]),1);

tmp = cellfun(@(X)size(X.model,2),corr_functions,'UniformOutput',false);

%- added a threshold of u
B = mat2cell(model.eigen_functions'*aggregate_mat*diag(double(u>0.5))*F,main_params.num_eigen,[tmp{:}]);



end

