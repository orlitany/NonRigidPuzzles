function show_eta_u(shape,parts,eta_u,eta_v,str)

N = numel(eta_u);
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:N
    subplot(2,N,i);showshape(shape,eta_u{i});caxis([0 1]);colorbar;colormap jet;title(str);    
    subplot(2,N,N+i);showshape(parts{i}.shape,eta_v{i});caxis([0 1]);colorbar;colormap jet;title(str);
end

end