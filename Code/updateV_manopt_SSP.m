function betta_hat = updateV_manopt_SSP(previous_betta,SSP_blk,C,Bu,sigma,parts,corr_functions,segment_areas,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


K = size(previous_betta{1},1);

previous_betta = vertcat(previous_betta{:});

for prt=1:main_params.num_parts - main_params.is_missing_part
    N = numel(parts{prt}.shape.X);
    tmp{prt} = sparse(diag2lin(C{prt}*parts{prt}.eigen_functions'*parts{prt}.M,corr_functions{prt}.parts));
    tmp{prt} = [tmp{prt} zeros(size(tmp{prt}))];
end
CA = blkdiag(tmp{:});

for prt=1:main_params.num_parts - main_params.is_missing_part
    N = numel(parts{prt}.shape.X);
    tmp{prt} = sparse(parts{prt}.shape.area_per_vertex');
    tmp{prt} = [tmp{prt} zeros(1,N)];
end
area_mat = blkdiag(tmp{:});

for prt=1:main_params.num_parts - main_params.is_missing_part
    tmp{prt} = sparse(gradientNorm(parts{prt}.shape,sigma{prt}));
    tmp{prt} = blkdiag(tmp{prt},tmp{prt});
end
regularizaion_mat = blkdiag(tmp{:});

for prt=1:main_params.num_parts - main_params.is_missing_part
    N = numel(parts{prt}.shape.X);
    tmp{prt} = repmat(speye(N),[1 2]);
end
aggregate_mat = blkdiag(tmp{:});

manifold = euclideanfactory(2*(main_params.num_parts-main_params.is_missing_part)*K);
problem.M = manifold;
% Phi = aggregate_mat'*model.eigen_functions;

Bu = horzcat(Bu{:});
Bu = Bu(:);
problem.cost  =  @(x)fun_v_SSP(x,SSP_blk,CA,Bu,regularizaion_mat,area_mat,aggregate_mat,segment_areas,algo_params.lambda_s_parts,algo_params.lambda_area_parts,algo_params.lambda_ones)
problem.egrad = @(x) fun_dv_SSP_fast(x,SSP_blk,CA,Bu,regularizaion_mat,area_mat,aggregate_mat,segment_areas,algo_params.lambda_s_parts,algo_params.lambda_area_parts,algo_params.lambda_ones)
                              

% Numerically check gradient consistency (optional).
% checkgradient(problem);
options.maxiter = algo_params.num_iter.manopt_maxiter;
% profile on 
betta  = conjugategradient(problem,previous_betta,options);
% profile viewer
betta_hat = reshape(betta,[K 2*(main_params.num_parts-main_params.is_missing_part)]);
betta_hat = num2cell(betta_hat,1);

end

