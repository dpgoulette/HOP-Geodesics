function cells1 = AlphaOneCells2d(DT,Edges,VV,VC)
%  AlphaOneCells2d - Calculates the alpha value for each Delaunay 1-cell.
%     The alpha that is returned for a particular edge (1-cell), is the
%     minumum value alpha must be to include that edge in the complex. This
%     alpha value comes from the study of alpha shapes.  Alpha shapes are
%     often used to reconstruct an object from a sample of points from that
%     object.  Alpha shapes (and the alpha value for each edge) allow us to
%     select a subset of the Delaunay 1-skeleton (the Delaunay graph) that
%     is based on the geometry of the data.
%
%  inputs:
%     DT - The Delaunay triangulation object for the data, DT.X.
%     Edges - The edges in the triangulation.
%     VV - The voronoi vertices for DT.X.
%     VC - The vertices of each voronoi cell surrounding each point in DT.X
%
%  output:
%     cells1 - k by 3 matrix.  A row in cells1 contains the indices a
%     Delaunay edge (in the first two columns) and the associated value of
%     alpha in the third column.  The index for an edge endpoint is an
%     index into the corresponding row in DT.X.

Data = DT.X;
cells1 = [Edges, zeros(length(Edges),1)];

for i=1:size(Edges,1)
    B=intersect(VC{Edges(i,1)}, VC{Edges(i,2)});
   
    if size(B,2) ~= 2
        fprintf('ERROR! This was unexpected.\n')
        s = strcat('Neighboring data-points in 2D must have',...
                   ' two Voronoi vertices in common. It is possible',...
                   ' that your data has a degenerate voronoi cell.');
        error(s)
    else
        % Get one of the data points from the current Delaunay edge.  (Both
        % endpoints are equidistant to their common voronoi wall so it doesn't
        % matter which one.)
        x=Data(Edges(i,1),:);
        % Get the two endpoints of the voronoi edge they share.
        e1=VV(B(1),:);
        e2=VV(B(2),:);
        
        v1=e2-e1;
        v2=x-e1;
        Theta1 = acos((v1*v2')/(norm(v1)*norm(v2)));
        
        w1=e1-e2;
        w2=x-e2;
        Theta2 = acos((w1*w2')/(norm(w1)*norm(w2)));
        
        % Determine if one of the endpoints of the voronoi edge is closest
        % to the current data points.  If not, then the closest point is
        % interior to the edge.
        if Theta1 > pi/2 %then e1 is closest point
            cells1(i,3) = norm(v2);
        elseif Theta2 > pi/2 %then e2 is closest point
            cells1(i,3) = norm(w2);
        else
            % The distance to the midpoint is closest.  So half the edge
            % length is alpha.  This means that the Delaunay edge between
            % the data points intersects the voronoi wall they share.
            cells1(i,3) = (1/2)*norm(Data(Edges(i,1),:)-Data(Edges(i,2),:));
        end
    end
end % main for loop
end % function
