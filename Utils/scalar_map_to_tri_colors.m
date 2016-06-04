function tri_colors = scalar_map_to_tri_colors(M, f, colorscale)

if all(f(1)==f)
    vertex_colors = colorscale(f,:);
else
    vertex_colors = colorscale(1+round((size(colorscale,1)-1)*(f-min(f))/range(f)),:);
end

m = size(M.TRIV,1);
tri_colors = zeros(m,3);
for i=1:m
    tri_colors(i,:) = mean(vertex_colors(M.TRIV(i,:),:));
end

end
