%%
clc, close all

load([main_params.output_folder '\' '/U step 5.mat']);

irr_prt = irrelevant_parts{1}.shape;
irr_prt.VERT = [irr_prt.X irr_prt.Y irr_prt.Z];

bd = calc_boundary_edges(irr_prt.TRIV);
A = irr_prt.VERT(bd(:,1),:);
B = irr_prt.VERT(bd(:,2),:);

figure, plot_mesh(irr_prt)

M = model.shape;
M.VERT = [M.X M.Y M.Z];
M.n = size(M.VERT,1);
M.m = size(M.TRIV,1);
M.S_tri = calc_tri_areas(M);

% unmatched points will receive white color
colors = [create_colormap(M,M); 1 1 1];

% figure
% plot_mesh(M), hold on
% plot_boundary_edges(M, calc_boundary_edges(M.TRIV))

M_seg = size(1+length(parts),1)*ones(M.n,1);

seg_colors = [...
    255 185 15 ; ...
    0 185 15 ; ...
    255 255 255] ./ 255;
%     0 185 15 ] ./ 255;
%     65 115 235 ; ...
%     255 185 15 ; ...
%     255 255 255]./255;

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
    
    
    
    bd = calc_boundary_edges(N.TRIV);
    A = N.VERT(bd(:,1),:);
    B = N.VERT(bd(:,2),:);
    A(A(:,3)>100,:)=[];
    B(B(:,3)>100,:)=[];
    
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
% 
% template_unmatched = true(M.n,1);
% for i=1:length(parts)
%     template_unmatched(eta_u{i}>0.5) = false;
% end
% M_seg(template_unmatched) = 3;
% 
% unmatched_part = (M_seg==3);
% corr_idx = (1:M.n)';
% corr_idx(unmatched_part) = size(colors,1);


%% show irrelevant part
N = irr_prt;
% N.VERT = [N.X N.Y N.Z];
N.n = size(N.VERT,1);
N.m = size(N.TRIV,1);
eta_v = zeros(N.n,1);
eta_u = zeros(M.n,1);
unmatched = eta_v<0.5;

% matches_masked = matches{i};
matches_masked(1:N.n) = size(colors,1);

% seg_masked = ones(N.n,1);
% seg_masked(unmatched) = 2;
% 
bd = calc_boundary_edges(N.TRIV);
A = N.VERT(bd(:,1),:);
B = N.VERT(bd(:,2),:);
A(A(:,3)>100,:)=[];
B(B(:,3)>100,:)=[];


figure
subplot(221), plot_scalar_map(M, eta_u), freeze_colors
subplot(222), plot_scalar_map(N, eta_v), freeze_colors
subplot(223), colormap(colors(1:end-1,:)), plot_scalar_map(M, 1:M.n), freeze_colors
subplot(224), colormap(colors(matches_masked,:)), plot_scalar_map(N, 1:N.n)

figure
plot_mesh(N), hold on
plot_boundary_edges(N, calc_boundary_edges(N.TRIV))
