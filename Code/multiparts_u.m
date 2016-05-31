function y = multiparts_u(fun, A, B, G, v, target_area, areas, mu1, mu2,lambda_ones, opt)
    
    num_parts = size(v,2);
    y = 0;
    for prt=1:num_parts - opt.is_missing_part
        
        y = y + fun(A{prt}, B, G{prt}, v(:,prt), ones(size(v,1),1), target_area{prt}, areas, mu1, mu2, opt);
                    
    end
    
    if opt.is_missing_part,
        y = y + fun(0*A{1},0*B,0*G{1},v(:,num_parts),ones(size(v,1),1),0,areas,0,mu2,opt);
    end
    
    y = y + 0.5*lambda_ones*sum( (sum(eta(v),2) - ones(size(v,1),1)).^2 );

end