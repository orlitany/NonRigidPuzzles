function M = part_to_M(model)

    M.TRIV = model.shape.TRIV;
    M.VERT = [model.shape.X model.shape.Y model.shape.Z];
    M.m = size(model.shape.TRIV,1);
    M.n = size(model.shape.X,1);
    M.S = model.M;  
    M.evecs = model.eigen_functions;
    
end