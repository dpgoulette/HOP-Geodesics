function cells1 = AlphaOneCells3d(DT,Edges,VV,VC)
%

cells1 = zeros(length(Edges),3);

% Loop through each edge in the Delaunay graph.
for i=1:size(Edges,1)
   % Find the voronoi vertices the two endpoints have in common.  These
   % voronoi vertices create a 2D planar polygon embedded in 3D.
   B=intersect(VC{Edges(i,1)}, VC{Edges(i,2)});
   Poly=VV(B,:);
   edge=DT.X(Edges(i,:),:);
   
   % Shift all points so that one edge enpoint is at the origin.
   Poly2=bsxfun(@minus,Poly,edge(1,:));
   e2=edge(2,:)-edge(1,:);
   
   % Use a Rodrigues rotation matrix to rotate the other end of the edge to
   % the x axis.  Thus the common polygonal wall will be parallel to the
   % y,z-axis.
   NormalVec = [0 e2(3) -e2(2)];        %cross product of e2 and [1 0 0]
   Theta=acos(e2(1)/norm(e2));          %angle between e2 and [1 0 0]
   RotMat = Rodrigues(Theta,NormalVec);
   PolyRot = RotMat*(Poly2)';           %rotate the polygon.
   
   PolyRot(1,:)=[];            %project onto the y-z plane. Poly is Kx2 now.
   PolyRot=PolyRot';
   
   % Now if the midpoint of the edge is already a vertex in the polygon,
   % then we are done with this edge.  The vertex is the midpoint of the
   % edge.  (In astronomical data with double floating point precision,
   % this case is practically impossible.)
   if ismember(PolyRot,[0 0],'rows') %midpoint was a vertex of Poly
      fprintf('Midpoint of edge %d was a vertex of its polygon!!!\n.',i)
      M= (edge(1,:)-edge(2,:))*(1/2);
      d= norm(edge(1,:)-M);
      cells1(i,:)=[Edges(i,:),d];
      % continue to the next edge.
      continue
   end
   
   
   % The key to the rest of this code: if the origin is interior to the
   % projected 2d polygon, then the midpoint of the original edge was
   % interior to the polygon.  If the origin is not interior then we need
   % to find the point on the boundary nearest to the origin.
   
   % Tack the origin on to the end of the polygon matrix.  Since the origin
   % is in the last row, the number of rows is equal to the origin row
   % index.
   MM=[PolyRot;0 0];
   originPolyID=size(MM,1);
   K=convhull(MM);
   
   % check to see if the origin is included in the convex hull.  If it is
   % not, then the origin is interior to the polygon.
   if ~(K == originPolyID)
      % then the midpoint of e1 and e2 was in the polygon to begin with
      % so we can get the point directly
      M= (edge(1,:)+edge(2,:))*(1/2); %midpoint
      d= norm(edge(1,:)-M); % alpha value is distance to midpoint
      cells1(i,:)=[Edges(i,:),d];
      
   elseif size(MM,1)==size(K,1)-1
      % Then the origin is not interior to the polygon AND the vertices in MM
      % are convex. No verts missed so we don't need to do new convhull
      % without the origin.
      
      O_index = find( K == originPolyID, 1, 'first' );
      
      % The points that precede and follow the origin in K are the closest
      % vertices to the origin. Note, because the origin is the last entry
      % in MM, it will never be first or last in K (because of the way
      % matlab does convhull).
      Vbefore=MM(K(O_index-1),:); %preceding vertex in K
      Vafter =MM(K(O_index+1),:); %following vertex in K
      % Get distances to the closest points.
      Nbefore = norm(Vbefore);
      Nafter = norm(Vafter);
      
      if Nbefore < Nafter
         V1 = Vbefore-Vafter;
         Angle = VecAngle(V1,Vbefore);
         if Angle > pi/2
            % Then Nbefore is the closest point to the origin.  So the
            % corresponding point in the original polygon in 3D was the
            % closest to the endpoints of the edge.
            
            d=norm(edge(1,:)-Poly(K(O_index-1),:));
            cells1(i,:)=[Edges(i,:),d]; %The most important line!
            
         else %Angle <= pi/2
            % The case when vbefore is closer with acute angle.  Since the
            % angle is not obtuse, the nearest point to the origin is in the
            % middle of the edge. Distance from Vbefor to the nearest point,
            % is L.  The variable t will hold the parameter needed to find
            % alpha.
            
            L = Nbefore*cos(Angle);
            t=L/norm(V1);
            
            % Find the real point in 3d corresponding to n in 2d. Then
            % find the distance from e1 to this point.  We do this with
            % a paramaterized line in 3D.
            n3= Poly(K(O_index-1),:) + (Poly(K(O_index+1),:) - ...
               Poly(K(O_index-1),:))*t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
            
         end
      else %Nbefore > Nafter
         % this is similar to the last case.
         V1 = Vafter-Vbefore;
         Angle = VecAngle(V1,Vafter);
         if Angle > pi/2
            d=norm(edge(1,:)-Poly(K(O_index+1),:));
            cells1(i,:)=[Edges(i,:),d];
            
         else %Angle <= pi/2
            %Angle is not obtuse so the nearest is in the
            %middle of the edge. Distance from Vafter to the nearest point,
            % is L.
            L = Nafter*cos(Angle);
            t=L/norm(V1);
            
            %find the real point in 3d corresponding to n in 2d. Then
            %find the distance from e1 to this point.  We do this with
            %a paramaterized line in 3D
            n3=Poly(K(O_index+1), :) + (Poly(K(O_index-1), :) -...
               Poly(K(O_index+1), :))*t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
         end
      end
      
   else
      % Then MM missed some of the original poly vertices vertices.  This
      % is the most complicated case.  There can now be multiple polygon
      % vertices with "direct line of sight" to the origin.  We need to
      % find the part of the polygon boundary that faces the origin.  It is
      % possible for a vertex on the opposite side of the polygon to be
      % closer to the origin than vertices that face the origin, but it
      % cannot be the closest point on the boundary (because the segement
      % connecting this vertex to the origin would cut through the center
      % of the polygon).
      
      % Find the index of the origin in K.
      O_index=find(K==originPolyID,1,'first');
      
      % Get the index of MM with preceding nbr of origin
      BeginningID = K(O_index-1);
      % shift the poly matrices (both 2d and 3d) so the point preceding the
      % origin is in the first row.  Matlab orders the convhull vertices
      % counterclockwise.
      Poly=circshift(Poly,-(BeginningID-1));
      PolyRot=circshift(PolyRot,-(BeginningID-1));
      
      Missed=size(MM,1)-size(K,1)+1;
      % Redo the convex hull without the origin.
      K2=convhull(PolyRot);
      
      % Get the indices of K that are on the boundary facing the origin.
      % Then get the points themselves.
      WallIDS=1:2+Missed;
      Wall=PolyRot(K2(WallIDS),:);
      
      % Distance to all facing wall vertices.
      WallNorms=sqrt(sum(Wall.^2,2));
      [WallMinD,WallMinID]=min(WallNorms);
      
      if WallMinID == 1
         V1=Wall(1,:)-Wall(2,:);
         Angle=VecAngle(V1,Wall(1,:));
         if Angle > pi/2
            % Then wall(1) is closest point on the boundary.
            d=norm(edge(1,:)-Poly(1,:));
            cells1(i,:)=[Edges(i,:),d];
         else
            % Nearest point is interior to the edge connecting wall(1) and
            % wall(2).
            L = WallMinD*cos(Angle);
            t = L/norm(V1);
            
            % Find the real point in 3d corresponding to n in 2d. Then
            % find the distance from e1 to this point.  We do this with
            % a paramaterized line in 3D
            n3=Poly(1,:)+(Poly(K2(2),:)-Poly(1,:))*t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
         end
         
      elseif WallMinID == length(WallIDS)
         % Similar to the last case but now the last vertex in the wall is
         % closest.
         V1=Wall(end,:)-Wall(end-1,:);
         Angle=VecAngle(V1,Wall(end,:));
         if Angle > pi/2
            % Then wall(end) is closest point on the boundary.
            
            d=norm(edge(1,:)-Poly(K2(WallIDS(end)),:));
            cells1(i,:)=[Edges(i,:),d];
         else
            L = WallMinD*cos(Angle);
            t = L/norm(V1);
            
            n3= Poly(K2(WallIDS(end)), :) + (Poly(K2(WallIDS(end-1)),:) - ...
                Poly(K2(WallIDS(end)), :)) * t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
         end
      else % The nearest wall vertex is in the middle of the wall
         V1=Wall(WallMinID,:)-Wall(WallMinID-1,:);
         V2=Wall(WallMinID,:)-Wall(WallMinID+1,:);
         Angle1=VecAngle(V1,Wall(WallMinID,:));
         Angle2=VecAngle(V2,Wall(WallMinID,:));
         
         if Angle1 > pi/2 && Angle2 > pi/2
            % Then Wall(WallMinID,:) is closest point on the boundary.
            
            d=norm(edge(1,:)-Poly(K2(WallIDS(WallMinID)),:));
            cells1(i,:)=[Edges(i,:),d];
         elseif Angle1 < pi/2
            % Then the nearest point is interior to the edge represented by
            % V1.
            L = WallMinD*cos(Angle1);
            t = L/norm(V1);
            
            n3 = Poly(K2(WallIDS(WallMinID)), :) + ...
                (Poly(K2(WallIDS(WallMinID-1)), :) - ...
                 Poly(K2(WallIDS(WallMinID)), :)) * t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
         else % Angele2 < pi/2
            % Then the nearest point is interior to the edge represented by
            % V2.
            L = WallMinD*cos(Angle2);
            t = L/norm(V2);
            
            n3= Poly(K2(WallIDS(WallMinID)), :) + ...
               (Poly(K2(WallIDS(WallMinID+1)), :) - ...
                Poly(K2(WallIDS(WallMinID)),:)) * t;
            d=norm(edge(1,:)-n3);
            cells1(i,:)=[Edges(i,:),d];
            
         end
      end
   end
   if mod(i,100) == 0
      fprintf('Finished %d out of %d edges.\n', i, size(Edges,1))
   end
end % main loop
end % main function

function RotMat = Rodrigues(Theta,u)
% Rodrigues - takes a vector axis of rotation and an angle of rotation. The
%       3-d rotation matrix is returned.  This is Rodrigues' rotation
%       formula.

u = u./norm(u,2);
S = [ 0    -u(3)  u(2);
      u(3)   0   -u(1);
     -u(2)  u(1)   0 ];
RotMat = eye(3) + sin(Theta)*S + (1-cos(Theta))*S^2;
end

function Theta = VecAngle(v1,v2)
%  Find the angle between two row vectors.

Theta = acos((v1*v2')/(norm(v1)*norm(v2)));
end
