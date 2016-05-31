function [X,Y,Z,TRIV,LT,v1,v2,vl] = edgecross(surface_x0, X,Y,Z,TRIV,LT,D,lt)

nv = length(X);
nt = size(TRIV,1);

t = find(lt(:,1) == lt(:,2) & lt(:,1) ~= lt(:,3));
tri = surface_x0.TRIV(t,:);
l1 = lt(t,1);
l2 = lt(t,3);

da = D(sub2ind(size(D), tri, repmat(l1,[1 3])));
db = D(sub2ind(size(D), tri, repmat(l2,[1 3])));

%lambda1 = (da(:,3)-da(:,1) + db(:,1) - db(:,3))./(db(:,1)-da(:,1));
%lambda2 = (da(:,3)-da(:,2) + db(:,2) - db(:,3))./(db(:,2)-da(:,2));

%delta_a = da(:,3)-da(:,1);
%delta_b = db(:,1)-db(:,3);
%lambda1 = (db(:,3) - da(:,1) + delta_b)./(delta_a + delta_b);
lambda1 = (db(:,1) - da(:,1)) ./ (da(:,3) - da(:,1) + db(:,1) - db(:,3));


%delta_a = da(:,3)-da(:,2);
%delta_b = db(:,2)-db(:,3);
%lambda2 = (db(:,3) - da(:,2) + delta_b)./(delta_a + delta_b);
lambda2 = (db(:,2) - da(:,2)) ./ (da(:,3) - da(:,2) + db(:,2) - db(:,3));

lambda1 = 1-lambda1;
lambda2 = 1-lambda2;


v1 = [surface_x0.X(tri(:,1)).*lambda1 + surface_x0.X(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Y(tri(:,1)).*lambda1 + surface_x0.Y(tri(:,3)).*(1-lambda1) , ...
      surface_x0.Z(tri(:,1)).*lambda1 + surface_x0.Z(tri(:,3)).*(1-lambda1)  ];

v2 = [surface_x0.X(tri(:,2)).*lambda2 + surface_x0.X(tri(:,3)).*(1-lambda2) , ...
      surface_x0.Y(tri(:,2)).*lambda2 + surface_x0.Y(tri(:,3)).*(1-lambda2) , ...
      surface_x0.Z(tri(:,2)).*lambda2 + surface_x0.Z(tri(:,3)).*(1-lambda2)  ];


X = [X; v1(:,1); v2(:,1)];
Y = [Y; v1(:,2); v2(:,2)];
Z = [Z; v1(:,3); v2(:,3)];
    
idx1 = [1:size(v1,1)] + nv;
idx2 = [1:size(v1,1)] + nv + size(v1,1);
idxt1 = [1:length(t)] + nt;
idxt2 = [1:length(t)] + nt + length(t);

TRIV(t,:) = [surface_x0.TRIV(t,1) surface_x0.TRIV(t,2) idx1(:)];
LT(t) = l1;
TRIV = [TRIV; [idx1(:) surface_x0.TRIV(t,2) idx2(:)]];
TRIV = [TRIV; [idx1(:) idx2(:) surface_x0.TRIV(t,3)]];
LT(idxt1) = l1;
LT(idxt2) = l2;

%vl = [t(:)];
vl = [l1(:) l2(:)];


