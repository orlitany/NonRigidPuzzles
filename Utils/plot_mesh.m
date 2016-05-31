function plot_mesh(N)
    trisurf(N.TRIV,N.VERT(:,1),N.VERT(:,2),N.VERT(:,3),zeros(size(N.VERT,1),1)),
    axis equal
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
    rotate3d on
end
