function maxclass = GeoTetras(maxclass)
%
%  ISSUE - GeoTetras uses the cell called "tetras" which we used to return.
%          This cell is no longer needed.  Consider removing and only using
%          maxclass.
%
% GeoTetras - searches through all of the geodesics in maxclass and finds
%     tetrahedra (in the graph theoretic sense).  This function is only
%     meaningful if the data is three dimensional.  A tetrahedron is
%     defined to be four maxima that are pairwise neighbors (i.e. there is
%     a geodesic connecting each pair) AND all six geodeiscs which make up
%     the "edges" of the tetrahedron are included (i.e. all six geodesics
%     were selected by the SelectGeodesics function).  A geodesic
%     tetrahedron will occur when four parwise adjacent maxima are
%     relatively close to each other (based on our geodesic metric).
%
%     input:
%        maxclass - the maxclass struct (with the geo_tris field empty).
%     output:
%        maxclass - the maxclass struct with the maxclass(i).geo_tetras field
%                   updated.
%
%     After GeoTetras completes, the maxclass.geo_tetras field will be
%     updated for each max.  If max i is not a member of any geodesic
%     tetrahedron, then maxclass(i).geo_tetras will be empty.  But if max i
%     is a member of at least one geodesic tetrahedron, then
%     maxclass(i).geo_tetras will contain a j by 3 cell (where j is the
%     number of tetrahedra that i is a member of).  Here are the contents
%     of the three columns:
%
%        maxclass(i).geo_tetras{:, 1}
%           contains the tetrahedron vertices.  These are the point indices of
%           the the four maxima that make up the tetrahedron (the index is
%           into DT.X)
%
%        maxclass(i).geo_tetras{:, 2}
%           contains the tetrahedron vertices.  These are the maxclass
%           indices of the the four maxima that make up the tetrahedron
%           (the index is into maxclass)
%
%        maxclass(i).geo_tetras{:, 3}
%           contains 1 or 0 depending if the tetrahedron is included or not
%           (respectively). A triangele is included if all six geodesics
%           making up the "sides" of the tetrahedron were included by the
%           SelectGeodesics function.  In this case we put a 1 in the
%           fourth column.  Otherwise we put a 0 to indicate the tetrahedron
%           is not included.

%%%%%%%%%%%%  MAIN FUNCTION %%%%%%%%%%%%%%

% We find each tetra attached to each max.  This is made efficient by
% starting our search with each triangle attached to each max.  A triangle
% plus one more vertex makes a tetra.

tetras = cell(length(maxclass),3);
for a = 1:length(maxclass)
   if isempty(maxclass(a).geo_tris)
      % So this max has no geodescis attached to it, hence it is not a part
      % of any tetra.
      continue
   end
   
   % loop through each triangle attached to max a
   for b = 1:size(maxclass(a).geo_tris, 1)
      
      %  Check to see of we already did this triangle.
      if maxclass(a).geo_tris{b, 2}(1) < a
         % this means we already cataloged all tetras attached to this
         % triangle in an earlier iteration.
         continue
      end
      
      current_tri = maxclass(a).geo_tris{b, 2};
      % Find any neighbors the 3 vertices have in common.  This defines a
      % tetra.
      nbs_1 = maxclass(current_tri(1)).nbormaxid;
      nbs_2 = maxclass(current_tri(2)).nbormaxid;
      nbs_3 = maxclass(current_tri(3)).nbormaxid;
      common_nbs = sort(intersect(intersect(nbs_1, nbs_2), nbs_3));
      
      % We know that all indices in current_tri are >= a. Now, if any of
      % the common neigbors have index less than a (the outer loop index),
      % then we already cataloged this tetra in an earlier iteration.  So
      % remove them from common_nbs.
      common_nbs(common_nbs < a) = [];
      
      % check if there are no common nbs
      if isempty(common_nbs)
         % Thus no common neighbors. So no tetras to catalog. Move to the
         % next.  Skip to the next triangle.
         continue
      end
      
      % Each common neighbor defines a tetra.  Store all tetras in the
      % Tetras cell for now.
      for c = 1:length(common_nbs)
         
         % we only want to catalog a tetra if it is in order.  This ensures
         % that we only catalog it one time.  We know current tri is in
         % order.  So only catalog if the current common neigbor index is
         % greater than all of the indices in current_tri.
         if common_nbs(c) < max(current_tri)
            continue
         end
         
         % Get the current tetra expressed both ways.
         % 1) with maxclass indices
         % 2) with raw data indices
         tet_maxclass_id = [current_tri, common_nbs(c)];
         tet_maxpoint_id = horzcat(maxclass(tet_maxclass_id).max);
         
         % Find out if the current tetra has all of its triangles included.
         %  This only occurs if all 6 of its geodesics are included.
         include_tet = tetra_membership(tet_maxclass_id, maxclass);
         
         % Each tetra is expressed with maxclass indices and also raw data
         % point indices.  Store these in the slot for each of the four
         % vertices.  (This work won't be repeated in later iterations.)
         tetras{a,2} = [tetras{a,2}; tet_maxclass_id];
         tetras{a,1} = [tetras{a,1}; tet_maxpoint_id];
         tetras{a,3} = [tetras{a,3}; include_tet];
         
         tetras{tet_maxclass_id(2), 2} = ...
            [tetras{tet_maxclass_id(2), 2}; tet_maxclass_id];
         tetras{tet_maxclass_id(2), 1} = ...
            [tetras{tet_maxclass_id(2), 1}; tet_maxpoint_id];
         tetras{tet_maxclass_id(2), 3} = ...
            [tetras{tet_maxclass_id(2), 3}; include_tet];
         
         tetras{tet_maxclass_id(3), 2} = ...
            [tetras{tet_maxclass_id(3), 2}; tet_maxclass_id];
         tetras{tet_maxclass_id(3), 1} = ...
            [tetras{tet_maxclass_id(3), 1}; tet_maxpoint_id];
         tetras{tet_maxclass_id(3), 3} = ...
            [tetras{tet_maxclass_id(3), 3}; include_tet];
         
         tetras{tet_maxclass_id(4), 2} = ...
            [tetras{tet_maxclass_id(4), 2}; tet_maxclass_id];
         tetras{tet_maxclass_id(4), 1} = ...
            [tetras{tet_maxclass_id(4), 1}; tet_maxpoint_id];
         tetras{tet_maxclass_id(4), 3} = ...
            [tetras{tet_maxclass_id(4), 3}; include_tet];
      end
   end
end

% Assign everything in tetras to the appropriate place in maxclass.
for a = 1:length(tetras)
   if ~isempty(tetras{a,1})
      maxclass(a).geo_tetras =...
         [num2cell(tetras{a, 1}, 2),...
         num2cell(tetras{a, 2}, 2),...
         num2cell(tetras{a, 3}, 2)];
   end
end
end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% DEPENDENT FUNCTIONS %%%%%%%%%%%%%%
%               | |       | |
%               | |       | |
%              \   /     \   /
%               \ /       \ /

function include_tet = tetra_membership(T, maxclass)
% tetra_membership -- checks to see if all 4 triangles of tetra T are
%        included.  This is the same as checking if all 6 geodesics that
%        make up T have been included.
%
%     inputs:   
%           T --- a 1 by 4 matrix of indices for the corners ofa tetra
%           maxclass --- The maxclass struct.
%
%     output:   
%           1 or 0 --- If the tetra is included return a 1. If not, return
%           a 0.

triangles = nchoosek(T, 3);
for a = 1:size(triangles,1)
   x = vertcat(maxclass(triangles(a,1)).geo_tris{:,2});
   r = getrow(x, triangles(a,:));
   % if we find that one of the triangles is not included, we are done.
   % Return 0.
   if ~maxclass(triangles(a,1)).geo_tris{r,4}
      include_tet = 0;
      % Short circut the function and return 0.
      return
   end
end

% If all of the triangles were included then we get to this point.
% Return 1
include_tet = 1;
end

function r = getrow(x, tri)
for i = 1:size(x,1)
   if isequal(x(i,:), tri)
      r = i;
      return
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


