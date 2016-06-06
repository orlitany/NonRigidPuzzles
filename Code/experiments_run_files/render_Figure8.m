%%
close all, clc

% params = get_faust_tx_params();
% params.resolution = [1600 1200];
% params.out_dir = './renderings_grey/';
% params.out_dir = './renderings/';

M = model.shape;
M.VERT = [M.X M.Y M.Z];
M.n = size(M.VERT,1);
M.m = size(M.TRIV,1);

Ns = cell(1,length(parts));
for i=1:length(parts)
    Ns{i} = parts{i}.shape;
    Ns{i}.VERT = [Ns{i}.X Ns{i}.Y Ns{i}.Z];
    Ns{i}.n = size(Ns{i}.VERT,1);
    Ns{i}.m = size(Ns{i}.TRIV,1);
end

color_res = 64;
colormaps = cell(1,length(parts));
% blue
colormaps{1} = [linspace(1,0,color_res)' linspace(1,0,color_res)' ones(color_res,1)];
% colormaps{1} = [ones(color_res,3)]; % this replaces the blue with white
%[255 100 0] is orange
colormaps{2} = [ones(color_res,1) linspace(1,100/255,color_res)' linspace(1,0/255,color_res)'];
% green
colormaps{3} = [linspace(1,0,color_res)' 1*ones(color_res,1) linspace(1,0,color_res)'];

% save_mesh_povray_with_normals(sprintf('%s/model.mesh',params.out_dir), M, ones(M.n,1), [1 1 1]);
% render_colors('model.mesh', params);
% copyfile(sprintf('%s/render_colors.png',params.out_dir), sprintf('%s/model.png',params.out_dir))

params.with_boundary = true;

for step=[1 2 8]
    
    load(sprintf('./Results/Figure8/U step %d.mat',step))
    load(sprintf('./Results/Figure8/V step %d.mat',step))
    
    figure
    subplot(231), plot_scalar_map(M, eta_u{1})
    subplot(232), plot_scalar_map(M, eta_u{2})
    subplot(233), plot_scalar_map(M, eta_u{3})
    subplot(234), plot_scalar_map(Ns{1}, eta_v{1})
    subplot(235), plot_scalar_map(Ns{2}, eta_v{2})
    subplot(236), plot_scalar_map(Ns{3}, eta_v{3})
    
end
%
colors = [create_colormap(M,M); 1 1 1];

%

for step=[1 2 8]    
    
    load(sprintf('./Results/Figure8/U step %d.mat',step))
    load(sprintf('./Results/Figure8/V step %d.mat',step))
    
    
    
    figure 
%     subplot(221), plot_scalar_map(M, eta_u{i}), freeze_colors
%     subplot(222), plot_scalar_map(Ns{i}, eta_v{i}), freeze_colors
    subplot(141), colormap(colors(1:end-1,:)), plot_scalar_map(M, 1:M.n), freeze_colors
    for i=1:3
    unmatched = eta_v{i}<0.5;
    matches_masked = matches{i};
    matches_masked(unmatched) = size(colors,1);
    subplot(1,4,i+1), colormap(colors(matches_masked,:)), plot_scalar_map(Ns{i}, 1:Ns{i}.n),freeze_colors
    end    
    
        
end

