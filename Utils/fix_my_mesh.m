function shape = fix_my_mesh(shape)
    M.TRIV = shape.TRIV;
    M.VERT = [shape.X shape.Y shape.Z];
    n=numel(shape.X);
    missing = setdiff(1:n, unique(shape.TRIV(:)));    
    is_outlier = false(n,1);
    is_outlier(missing) = true;
    if sum(is_outlier)>0
        M = removeVertices(M, is_outlier, false);
    end
    


end