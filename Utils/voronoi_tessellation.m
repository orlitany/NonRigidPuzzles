function [surface_x1,edges,D] = voronoi_tessellation(surface_x0, D)

    [~,l] = min(D,[],2);
    
    
lt = l(surface_x0.TRIV);

t = find(lt(:,2) == lt(:,3) & lt(:,1) ~= lt(:,2));
surface_x0.TRIV(t,:) = surface_x0.TRIV(t,[2 3 1]);

lt = l(surface_x0.TRIV);
t = find(lt(:,1) == lt(:,3) & lt(:,1) ~= lt(:,2));
surface_x0.TRIV(t,:) = surface_x0.TRIV(t,[3 1 2]);

lt = l(surface_x0.TRIV);

TRIV = surface_x0.TRIV;
X = surface_x0.X;
Y = surface_x0.Y;
Z = surface_x0.Z;
LT = median(l(surface_x0.TRIV),2);

[X,Y,Z,TRIV,LT,v1,v2,vl] = edgecross(surface_x0, X,Y,Z,TRIV,LT,D,lt);
%[X,Y,Z,TRIV,LT,v1_,v2_,v3_,v,vl_] = vertcross(surface_x0, X,Y,Z,TRIV,LT,D,lt);
[X,Y,Z,TRIV,LT,v1_,v2_,v3_,v,el] = vertcross(surface_x0, X,Y,Z,TRIV,LT,D,lt);

surface_x1 = surface_x0;
surface_x1.TRIV = TRIV;
surface_x1.X = X;
surface_x1.Y = Y;
surface_x1.Z = Z;
surface_x1.tri_labels = LT;

edges = [];
edges.X = [v1(:,1)'; v2(:,1)'];
edges.Y = [v1(:,2)'; v2(:,2)'];
edges.Z = [v1(:,3)'; v2(:,3)'];
edges.label = [vl];

edges.X = [edges.X [v1_(:,1)'; v(:,1)']];
edges.X = [edges.X [v2_(:,1)'; v(:,1)']];
edges.X = [edges.X [v3_(:,1)'; v(:,1)']];
%edges.label = [edges.label; vl_; vl_; vl_;];
edges.label = [edges.label; el;];

edges.Y = [edges.Y [v1_(:,2)'; v(:,2)']];
edges.Y = [edges.Y [v2_(:,2)'; v(:,2)']];
edges.Y = [edges.Y [v3_(:,2)'; v(:,2)']];

edges.Z = [edges.Z [v1_(:,3)'; v(:,3)']];
edges.Z = [edges.Z [v2_(:,3)'; v(:,3)']];
edges.Z = [edges.Z [v3_(:,3)'; v(:,3)']];
