function bd = calc_boundary_edges(triangles)

if isstruct(triangles)
    triangles = triangles.TRIV;
end

[c,d,~,~] = get_boundary(triangles);
bd = zeros(length(d),2);

for i=1:length(c)
    t = triangles(c(i),:);
    v = true(1,3);
    v(d(i)) = false;
    v = t(v);
    bd(i,1) = v(1);
    bd(i,2) = v(2);
end

function [ c,d,I,v ] = get_boundary( tri )
%GET_BOUNDARY determines the boundary edges of a triangular mesh 
%   [c,d] = get_boundary(tri) takes as input a list tri of consistently oriented triangles
%   returns the indices c of the triangles the boundary edges belong to and
%   the (local) indices d (in {1,2,3}) of the vertices opposing the boundary
%   edge
%   One gets the global indices of those vertices via F(sub2ind(size(F),c,d))
%   Via 
%   d1 = mod(d+1,3); d1(d1==0) = 3;
%   d2 = mod(d+2,3); d2(d2==0) = 3;
%   one gets the local indices of the boundary vertices.


% Transpose matrix if neccessary
if size(tri,1)<size(tri,2)
    tri=tri';
end
m = size(tri,1);

%% Check for opposing halfedges

% Matrix of directed edges
I = [tri(:,1) tri(:,2);
     tri(:,2) tri(:,3);
     tri(:,3) tri(:,1)];

b = not(ismember(I(:,[2 1]),I,'rows'));
b = find(b);



% Triangle indices
c = mod(b,m);
c(c==0) = m;

% vertex opposing boundary edge
d = floor((b-1)/m);
d(d==0)=3;


% % Directed boundary edges
I=I(b,:);

% Boundary vertices
[~,~,v] = find(I);
v = unique(v);

end

end
