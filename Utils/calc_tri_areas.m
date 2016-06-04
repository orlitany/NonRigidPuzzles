function S_tri = calc_tri_areas(M)

S_tri = zeros(size(M.TRIV,1),1);

for k=1:size(M.TRIV,1)
    e1 = M.VERT(M.TRIV(k,3),:) - M.VERT(M.TRIV(k,1),:);
    e2 = M.VERT(M.TRIV(k,2),:) - M.VERT(M.TRIV(k,1),:);
    S_tri(k) = 0.5*norm(cross(e1,e2));
end

end
