function [v_hat,A] = updateV(C,B,sigma,parts,corr_functions,segment_areas,main_params,algo_params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for prt = 1:main_params.num_parts - main_params.is_missing_part
    
    N = numel(parts{prt}.shape.X);
    G_part = corr_functions{prt}.parts;
    regularizaion_mat = gradientNorm(parts{prt}.shape,sigma{prt});
    

    cvx_begin
    
        variable v(N) nonnegative
        minimize( 0.5*square_pos( norm( C{prt}*parts{prt}.eigen_functions'*parts{prt}.M*diag(v)*G_part - B{prt} ,'fro')) +...
                  0.5*algo_params.lambda_s_parts*v'*regularizaion_mat*v + ... 
                 0.5*algo_params.lambda_area*sum( parts{prt}.shape.area_per_vertex'*v - segment_areas(prt) ).^2 ) 
    
        subject to
            v <= 1
         
    cvx_end
    figure;showshape(parts{prt}.shape,v);title('v');caxis([0 1]);colorbar
    
    v_hat{prt} = v;
    A{prt} = parts{prt}.eigen_functions'*parts{prt}.M*diag(v)*G_part;
    
    
end





end

