%% Add paths and prepare settings

%%
fprintf('Loading shapes...');tStart = tic;
S1 = load('./data/mesh000.mat');S1.name = 'mesh000.mat';
S2 = load('./data/mesh002.mat');S2.name = 'mesh002.mat';
fprintf('done, took %.1f \n',toc(tStart))
%%
do_plot = true;
save_pth = './data_resized/';

for ii = 1:4
    nVertTgt = ii * 1e3;
    remes_opts = struct('placement',0,'vertices',nVertTgt,'verbose',0);
    [TRIV_, X_, Y_, Z_] = remesh(S1.shape, remes_opts);
    
    % find original indices
    verts = [S1.shape.X S1.shape.Y S1.shape.Z ];
    dist_from_centroid = sqrt(sum(bsxfun(@minus,verts,mean(verts)).^2,2));
    [orig_idx,d] = knnsearch(verts,[X_, Y_, Z_]);
    
    assert (max(d)/max(dist_from_centroid) < 1e5,'distanced seem high, perhaps remesh failed')
    
    
    cur_save_pth = [save_pth int2str(ii) 'K/'];
    MakeDir(cur_save_pth)
    
    surface = struct('X',X_,'Y',Y_,'Z',Z_,'TRIV',TRIV_,'orig_idx',orig_idx); %#ok<NASGU>
    save([cur_save_pth S1.name],'surface')
    if do_plot
        figure(2*ii);clf;trisurf(surface.TRIV,surface.X,surface.Y,surface.Z);axis equal
    end
    
    surface = struct(...
        'X',S2.shape.X(orig_idx),...
        'Y',S2.shape.Y(orig_idx),...
        'Z',S2.shape.Z(orig_idx),'TRIV',TRIV_,'orig_idx',orig_idx);
    save([cur_save_pth S2.name],'surface')
    if do_plot        
        figure(2*ii-1);clf;trisurf(surface.TRIV,surface.X,surface.Y,surface.Z);axis equal
    end
    
end


