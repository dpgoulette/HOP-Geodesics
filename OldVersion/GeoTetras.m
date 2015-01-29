function [tetras, maxclass] = GeoTetras(maxclass)
%

%%%%%%%%%%%%  MAIN FUNCTION %%%%%%%%%%%%%%
tetras = cell(length(maxclass),3);
for a = 1:length(maxclass)
   if isempty(maxclass(a).geo_tris)
      % So this max is not a part of any max triangle.  Thus it can't be a
      % part of a tetra. Move on to the next max.
      continue
   end
   for b = 1:size(maxclass(a).geo_tris, 1)
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
      % then we already cataloged this tetra in an earlier iteration.
      common_nbs(common_nbs < a) = [];
      
      % check if there are no common nbs
      if isempty(common_nbs)
         % Thus no common neighbors. So no tetras to catalog. Move to the
         % next.
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

% % Sort the tetras?? Do we need this?
% for a = 1:length(tetras)
%    if ~isempty(tetras{a,1})
%       tetras{a,1} = sortrows(tetras{a,1});
%       tetras{a,2} = sortrows(tetras{a,2});
%    end
% end
for a = 1:length(tetras)
   if ~isempty(tetras{a,1})
      maxclass(a).geo_tetras =...
         [num2cell(tetras{a, 1}, 2),...
         num2cell(tetras{a, 2}, 2),...
         num2cell(tetras{a, 3}, 2)];
   end
end
end % main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% Dependent functions for GeoTetras %%%%%
function include_tet = tetra_membership(T, maxclass)
% this function checks to see if all 4 triangles of tetra T are
% included.  This is the same as checking if all 6 geodesics that
% make up T have been included.
%
% inputs:   T        --- 1 by 4 matrix of indices for the corners of
%                        a tetra
%           maxclass --- The maxclass struct.
%
% output:   1 or 0   --- If the tetra is included return a 1.
%                        If not, return a 0

% example T = [1, 6, 10, 51]
triangles = nchoosek(T, 3);
for a = 1:size(triangles,1)
   x = vertcat(maxclass(triangles(a,1)).geo_tris{:,2});
   r = getrow(x, triangles(a,:));
   % if we find that one of the triangles is not included, we are done.
   % Return 0.
   if ~maxclass(triangles(a,1)).geo_tris{r,4}
      include_tet = 0;
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


