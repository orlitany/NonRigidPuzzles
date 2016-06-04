function save_mesh_povray_with_normals(fname, M, colors, colorscale)

if size(colors,1)==1
    colors = colors';
end

if size(colors,1) == M.n
    colors = scalar_map_to_tri_colors(M, colors, colorscale);
elseif size(colors,1) ~= M.m
    error('The given scalar function should be defined either on vertices or triangles.');
end

normals = calc_normals(M.VERT', M.TRIV')';

fid = fopen(fname, 'w');

n = size(M.VERT,1);
m = size(M.TRIV,1);

fprintf(fid, '%d,%d,', n, m);

X = M.VERT;
X(isnan(X(:,1)),:) = 0;
fprintf(fid, '%.4f,%.4f,%.4f,', X');

X = normals;
X(isnan(X(:,1)),:) = 0;
fprintf(fid, '%.4f,%.4f,%.4f,', X');

fprintf(fid, '%d,%d,%d,', (M.TRIV-1)');

fprintf(fid, '%.4f,%.4f,%.4f,', colors');

fclose(fid);

end
