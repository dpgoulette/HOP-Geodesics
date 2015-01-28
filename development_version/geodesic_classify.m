function maxclass = geodesic_classify(maxclass)
% geodesic_classify - adds a 5th column to the cell in maxclass(i).geodesics.  This column will
%     This column will hold a 0, 1, 2 or 3.
%        0 if the geodesic is not included
%        1 if the geodesic is included but not a member of any triangle
%        2 if the geodesic is included in a triangle but not in any tetra
%        3 if teh geodesic is included in a tetra (only applicable in 3D)
%
%     (Note that the maxclass.geo_tetras field will be empty if the data is
%     2D so there won't be a 3 for any geodesic in that case.)  These
%     numbers correspond to the "demension" of the local space in an
%     approximate sense.  If a tetradedron has a 1, then it is in a roughly
%     1-dimensional structure.  If it has a 2 then it is in a roughly two
%     dimensional structure (a triangle embedded in either 2D or 3D space).
%     If it has a 3 then it is roughly in a 3D local structure embedded in
%     3D.
%
%     input:
%        maxclass - the maxclass struct.
%     output:
%        maxclass - the maxclass struct with the fifth column of
%                   maxclass.geodesics added.


for a = 1: length(maxclass)
   
   % If there are no geodesics connected to the max there is nothing to do.
   % So continue to the next max.
   if isempty(maxclass(a).geodesics);
      continue
   end
   
   % Get the geodesic endpoints connected to the current max
   geos = vertcat(maxclass(a).geodesics{:,2});
   
   % Assign NaN to the third column entries. We will use this slot for
   % tracking progress.  Mark zero in the third column for any geodesic not
   % included.
   geos(:,3) = vertcat(maxclass(a).geodesics{:,4});
   geos(geos(:,3) == 2, 3) = NaN;
   geos(geos(:,3) == 1, 3) = 0;
   
   % Now check if this max is a member of any tetras.  
   if ~isempty(maxclass(a).geo_tetras) &&...
         any(vertcat(maxclass(a).geo_tetras{:,3}));
      % Then there are tetras.  And at least one of them is included.
      % Find all of the geodesics which are included in tetras, and include
      % the current max point.  They need to be expressed in the order that
      % they are in maxclass(i).geodesics{j,2}(:,[1,2])
      
      % Find which geodesics are a member of a tetra.
      geos = find_tet_geos(geos, maxclass, a);
      
      % Find which of the remaining geodesics are members of a triangle
      % (the remaining geos will have a NaN).
      if any(isnan(geos(:,3)))
         geos = find_tri_geos(geos, maxclass, a);
      end
   elseif ~isempty(maxclass(a).geo_tris) &&...
         any(vertcat(maxclass(a).geo_tris{:,4}));
      % Thus there were no tetras connected to this max but there are
      % triangles included.
      geos = find_tri_geos(geos, maxclass, a);
   end
   
   % Any leftover NaNs must be geodesics that are included but not a member
   % of any tetra or triangle.  Put a 1.
   if any(isnan(geos(:,3)))
      geos(isnan(geos)) = 1;
   end
   
   % Store the results in maxclass(a).geodesics
   maxclass(a).geodesics(:,5) = num2cell(geos(:,3));
   
end
end % main function: geodesic_classify

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% DEPENDENT FUNCTIONS %%%%%%%%%%%%%%
%               | |       | |
%               | |       | |
%              \   /     \   /
%               \ /       \ /

function geos = find_tet_geos(geos, maxclass, a)
% find_tet_geos -- finds any geodesics connected to max a which are members
%     of a geodesic tetrahedron.  The third column of the geos matrix is
%     updated.  A 3 is entered if the geodesic in that row is a member
%     of some tetra.
%
%     inputs:
%              a           - the index of a max
%              maxclass
%              geos        - the geos matrix (see parent function)
%
%     output:
%              geos        - updated with a 3 in the third column for each
%                            geodesic that is a member of a tetra.


% get the current max's tetras that are included
temp1 = vertcat(maxclass(a).geo_tetras{:,3});
in = temp1 == 1;
temp = vertcat(maxclass(a).geo_tetras{in,1});
% get the unique vertices
temp = unique(temp);
% remove the current max.  we want the neighbors
temp(temp == maxclass(a).max) = [];
% find the geos that are in tetras
r = ismember(geos(:,2), temp);
% assign 3 to those rows.  They are 3D geos
geos(r,3) = 3;
end % find_tet_geos

function geos = find_tri_geos(geos, maxclass, a)
% find_tri_geos -- finds any geodesics connected to max a which are members
%     of a geodesic triangle.  The third column of the geos matrix is
%     updated.  A 2 is entered if the geodesic in that row is a member
%     of some triangle.
%
%     inputs:
%              a           - the index of a max
%              maxclass
%              geos        - the geos matrix (see parent function)
%
%     output:
%              geos        - updated with a 2 in the third column for each
%                            geodesic that is a member of a triangle.


% get the current max's tris
temp1 = vertcat(maxclass(a).geo_tris{:,4});
in = temp1 == 1;
temp = vertcat(maxclass(a).geo_tris{in,1});
% get the unique vertices
temp = unique(temp);
% remove the current max.  we want the neighbors
temp(temp == maxclass(a).max) = [];
% find the geos that are in triangles
r = ismember(geos(:,2), temp) & isnan(geos(:,3));
% assign 2 to those rows.  They are 2D geos
geos(r,3) = 2;
end % find_tri_geos