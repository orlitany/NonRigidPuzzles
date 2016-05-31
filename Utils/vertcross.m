function [X,Y,Z,TRIV,LT,v1,v2,v3,v4,vl] = vertcross(surface_x0, X,Y,Z,TRIV,LT,D,lt)

nv = length(X);
nt = size(TRIV,1);

t = find(lt(:,1) ~= lt(:,2) & lt(:,1) ~= lt(:,3) & lt(:,2) ~= lt(:,3));
tri = surface_x0.TRIV(t,:);
l1 = lt(t,1);
l2 = lt(t,2);
l3 = lt(t,3);

da = D(sub2ind(size(D), tri, repmat(l1,[1 3])));
db = D(sub2ind(size(D), tri, repmat(l2,[1 3])));
dc = D(sub2ind(size(D), tri, repmat(l3,[1 3])));

lambda1 = 1-(dc(:,1) - da(:,1)) ./ (da(:,3) - da(:,1) + dc(:,1) - dc(:,3));
lambda2 = 1-(db(:,1) - da(:,1)) ./ (da(:,2) - da(:,1) + db(:,1) - db(:,2));
lambda3 = 1-(dc(:,2) - db(:,2)) ./ (db(:,3) - db(:,2) + dc(:,2) - dc(:,3));

v1 = [surface_x0.X(tri(:,1)).*lambda1 + surface_x0.X(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Y(tri(:,1)).*lambda1 + surface_x0.Y(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Z(tri(:,1)).*lambda1 + surface_x0.Z(tri(:,3)).*(1-lambda1)  ];

v2 = [surface_x0.X(tri(:,1)).*lambda2 + surface_x0.X(tri(:,2)).*(1-lambda2) , ...
      surface_x0.Y(tri(:,1)).*lambda2 + surface_x0.Y(tri(:,2)).*(1-lambda2) , ...
      surface_x0.Z(tri(:,1)).*lambda2 + surface_x0.Z(tri(:,2)).*(1-lambda2)  ];

v3 = [surface_x0.X(tri(:,2)).*lambda3 + surface_x0.X(tri(:,3)).*(1-lambda3) , ...
      surface_x0.Y(tri(:,2)).*lambda3 + surface_x0.Y(tri(:,3)).*(1-lambda3) , ...
      surface_x0.Z(tri(:,2)).*lambda3 + surface_x0.Z(tri(:,3)).*(1-lambda3)  ];
  
v4 = (v1+v2+v3)/3;  
  
X = [X; v1(:,1); v2(:,1); v3(:,1); v4(:,1)];
Y = [Y; v1(:,2); v2(:,2); v3(:,2); v4(:,2)];
Z = [Z; v1(:,3); v2(:,3); v3(:,3); v4(:,3)];

    
idx1 = [1:size(v4,1)] + nv;
idx2 = [1:size(v4,1)] + nv + size(v4,1);
idx3 = [1:size(v4,1)] + nv + size(v4,1)*2;
idx4 = [1:size(v4,1)] + nv + size(v4,1)*3;


idxt2 = [1:length(t)] + nt;
idxt3 = [1:length(t)] + nt + length(t);
idxt4 = [1:length(t)] + nt + length(t)*2;
idxt5 = [1:length(t)] + nt + length(t)*3;
idxt6 = [1:length(t)] + nt + length(t)*4;

TRIV(t,:) = [surface_x0.TRIV(t,1) idx2(:) idx4(:)];     % 1  
LT(t) = l1;
TRIV = [TRIV; [idx2(:) surface_x0.TRIV(t,2) idx4(:)]];  %2
TRIV = [TRIV; [surface_x0.TRIV(t,1) idx4(:) idx1(:)]];  %3 
TRIV = [TRIV; [surface_x0.TRIV(t,2) idx3(:) idx4(:)]];  %4
TRIV = [TRIV; [idx1(:) idx4(:) surface_x0.TRIV(t,3)]];  %5
TRIV = [TRIV; [idx4(:) idx3(:) surface_x0.TRIV(t,3)]];  %6

LT(idxt2) = l2;
LT(idxt3) = l1;
LT(idxt4) = l2;
LT(idxt5) = l3;
LT(idxt6) = l3;

%vl = [t(:)];
%vl = [l1(:) l2(:) l3(:)];
vl = [ [l1(:) l3(:)];  [l1(:) l2(:)];  [l2(:) l3(:)] ];

