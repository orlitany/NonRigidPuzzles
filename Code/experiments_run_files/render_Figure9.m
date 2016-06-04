%%
clc, close all

load([main_params.output_folder '\' 'V step 8.mat']);

M = model.shape;
M.VERT = [M.X M.Y M.Z];
M.n = size(M.VERT,1);
M.m = size(M.TRIV,1);
[~,~,~,M.S_tri] = extract_eigen_functions_new(M,1);

% D = calc_dist_matrix(M);

% unmatched points will receive white color
colors = [create_colormap(M,M); 1 1 1];

figure
plot_mesh(M), hold on
plot_boundary_edges(M, calc_boundary_edges(M.TRIV))

M_seg = size(1+length(parts),1)*ones(M.n,1);

seg_colors = [...
    165 226 117 ; ...
    255 160 255;...
    0 185 15 ; ...
    65 115 235 ; ...
    255 185 15 ; ...
    255 255 255]./255;

errors = cell(1,length(parts));

for i=1:length(parts)
    
    N = parts{i}.shape;
    N.VERT = [N.X N.Y N.Z];
    N.n = size(N.VERT,1);
    N.m = size(N.TRIV,1);
    
    % make sure there are matches for all points
    assert(N.n == length(matches{i}));
    
    unmatched = eta_v{i}<0.5;
    
    matches_masked = matches{i};
    matches_masked(unmatched) = size(colors,1);
    
    seg_masked = ones(N.n,1);
    seg_masked(unmatched) = 2;
    
    fname = sprintf('part%d',i);
    
    
    M_seg(eta_u{i}>0.5) = i;
    
    figure
    subplot(221), plot_scalar_map(M, eta_u{i}), freeze_colors
    subplot(222), plot_scalar_map(N, eta_v{i}), freeze_colors
    subplot(223), colormap(colors(1:end-1,:)), plot_scalar_map(M, 1:M.n), freeze_colors
    subplot(224), colormap(colors(matches_masked,:)), plot_scalar_map(N, 1:N.n)
    
    figure
    plot_mesh(N), hold on
    plot_boundary_edges(N, calc_boundary_edges(N.TRIV))
    
end

%%