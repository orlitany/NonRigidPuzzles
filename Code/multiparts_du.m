function dy = multiparts_du(fun, A, B, G, v, target_area, areas, mu1, mu2,lambda_ones, opt)
    
   num_parts = size(v,2);
   
   dy = zeros(size(v,1),num_parts);
   sum_v = lambda_ones*(sum(eta(v),2) - ones(size(v,1),1));
   
   for prt=1:num_parts - opt.is_missing_part
       dy(:,prt) = fun(A{prt}, B, G{prt}, v(:,prt), ones(size(v,1),1), target_area{prt}, areas, mu1, mu2, opt,1) + ...
                    sum_v.*diff_eta(v(:,prt));
       
   end
   if opt.is_missing_part,
       dy(:,num_parts) = fun(A{1},B,G{1}, v(:,num_parts), ones(size(v,1),1), 0, areas, 0, mu2, opt,0) + ...
           sum_v.*diff_eta(v(:,num_parts));
   end
   
   
end