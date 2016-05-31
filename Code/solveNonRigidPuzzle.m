function [C,eta_u,eta_v] = solveNonRigidPuzzle(model,parts,corr_functions,main_params,algo_params)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
rng(1)

%- initialize
num_points = numel(model.shape.X);

v = cellfun(@(P)ones(numel(P.shape.X),1),parts,'UniformOutput',false);
eta_v = cellfun(@(P)ones(numel(P.shape.X),1),parts,'UniformOutput',false);

[model.eigen_functions,~,model.eigen_values,~] = extract_eigen_functions_new(model.shape,main_params.num_eigen);

[~,~,~,model.M] = extract_eigen_functions_new(model.shape,10);
model.shape.area_per_vertex = full(diag(model.M));

for prt = 1:main_params.num_parts        
    [parts{prt}.eigen_functions,~, parts{prt}.eigen_values,parts{prt}.M] = extract_eigen_functions_new(parts{prt}.shape,main_params.num_eigen);    
    parts{prt}.shape.area_per_vertex = full(diag(parts{prt}.M));
    A{prt} = parts{prt}.eigen_functions'*parts{prt}.M*(diag(v{prt}) * corr_functions{prt}.parts);
    u{prt} = ones(num_points,1);
    B{prt} = model.eigen_functions'*model.M*(diag(u{prt}) * corr_functions{prt}.model);    
    est_rank{prt} = sum(parts{prt}.eigen_values - max(model.eigen_values)<0);        
    W{prt} = calcW(main_params.num_eigen,est_rank{prt});    
    d{prt} = zeros(1,main_params.num_eigen); d{prt}(1:est_rank{prt})=1;    
    D{prt} = repmat(d{prt},main_params.num_eigen,1);
end

M = part_to_M(model);
N = cellfun(@(X)part_to_M(X),parts,'UniformOutput',false);

%- initialize the parts areas as all 1's.
for prt = 1:main_params.num_parts
    parts_areas{prt} = sum(parts{prt}.shape.area_per_vertex);%model.shape.area_per_vertex'*double(model.shape.ind_labels==prt); 
end
opt.is_missing_part = main_params.is_missing_part;

eta_u_null = ones(num_points,1);

C = cellfun(@(X)(max(max(X))-X)./max(max(X)),W,'UniformOutput',false);
for C_iter = 1:algo_params.num_iter.C_step + algo_params.num_iter.ransac   
    %- C-step (and icp)
    C = updateC21(C,A,B,W,d,D,main_params,algo_params);           
    parfor prt=1:main_params.num_parts
        [Co{prt},matches{prt}] = run_icp(M, N{prt}, est_rank{prt}, C{prt}, algo_params.num_iter.icp)
    end    
    C = cellfun(@(X,Y)[X zeros(size(Y,1),size(Y,1)-size(X,2))],Co,C,'UniformOutput',false);        
    
    if (C_iter > algo_params.num_iter.C_step)
        for  prt=1:main_params.num_parts
            idx=fps_euclidean([parts{prt}.shape.X parts{prt}.shape.Y parts{prt}.shape.Z], algo_params.num_ransac_corr, randi(numel(parts{prt}.shape.X)));
            FG = compute_indicator_functions({M,N{prt}}, [matches{prt}(idx), idx]', 0.1);
            corr_functions{prt}.model = FG{1};
            corr_functions{prt}.parts = FG{2};   
            A{prt} = parts{prt}.eigen_functions'*parts{prt}.M*(diag(v{prt}) * corr_functions{prt}.parts);
        end        
    end
    
    %- U-step    
    eta_u = cellfun(@(X,P,Y)model.eigen_functions*X*(P.eigen_functions'*(P.M*Y)),C,parts,v,'UniformOutput',false);
    u = eta_u;
    if main_params.is_missing_part         
        u{main_params.num_parts+1} = double(sum(cell2mat(eta_u),2)<0.1);%eta_u_null;                 
    end 
       
        
    opt.tv_sigma = 0.2; 
    opt.tv_mean = 0.5;      
    u = optimize_u_multi_parts(u,M,N,corr_functions,C,parts_areas,algo_params.lambda_area_model,algo_params.lambda_reg_model,algo_params.lambda_ones,opt,algo_params.num_iter.manopt_maxiter_u);
    eta_u = cellfun(@(X)eta(X),u(1:main_params.num_parts),'UniformOutput',false);
    if main_params.is_missing_part,
        eta_u_null = eta(u{main_params.num_parts+1});
    end
    
    save([main_params.output_folder 'U step ' num2str(C_iter)],'C','eta_u','eta_v','matches');
    
    B = cellfun(@(x,F)model.eigen_functions'*(model.M*(diag(x)*F.model)),eta_u,corr_functions,'UniformOutput',false);
    
    %- V-step
    if main_params.is_clutter,
    
        parfor prt=1:main_params.num_parts
            v{prt} = optimize_v_parts(M,N{prt},corr_functions{prt},C{prt},eta_u{prt},algo_params.lambda_area_parts,algo_params.lambda_reg_parts,opt,algo_params.num_iter.manopt_maxiter_v);        
            eta_v{prt} = eta(v{prt});
            parts_areas{prt} = parts{prt}.shape.area_per_vertex'*eta_v{prt};
            A{prt} = parts{prt}.eigen_functions'*(parts{prt}.M*(diag(v{prt}) * corr_functions{prt}.parts));
        end
        save([main_params.output_folder 'V step ' num2str(C_iter)],'C','eta_u','eta_v','matches');
    end
    
    show_eta_u(model.shape,parts,eta_u,eta_v,['step number ' num2str(C_iter)]);drawnow();

end
    
end

