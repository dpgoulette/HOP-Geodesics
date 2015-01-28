function maxclass = GeoTris(maxclass)
%
%  ISSUE - This code uses a cell Geodesic_Tris which we no longer use (and
%          the function doesn't return).  We now only use maxclass.
%          Consider removing Geodesic_Tris from this code.
%
% GeoTris - searches through all of the geodesics in maxclass and finds
%     triangles (in the graph theoretic sense).  A triangle is defined to
%     be three maxima that are pairwise neighbors (i.e. there is a geodesic
%     connecting each pair) AND all three geodeiscs are included (i.e. all
%     three were selected by the SelectGeodesics function).  A geodesic
%     triangle will occur when three parwise adjacent maxima are relatively
%     close to each other (based on our geodesic metric).
%
%     input:
%        maxclass - the maxclass struct (with the geo_tris field empty).
%     output:
%        maxclass - the maxclass struct with the maxclass(i).geo_tris field
%                   updated.
%
%     After GeoTris completes, the maxclass.geo_tris field will be updated
%     for each max.  If max i is not a member of any geodesic triangle,
%     then maxclass(i).geo_tris will be empty.  But if max i is a member of
%     at least one geodesic triangle, then maxclass(i).geo_tris will
%     contain a j by 4 cell (where j is the number of triangles that i is a
%     member of).  Here are the contents of the four columns:
%
%        maxclass(i).geo_tris{:, 1}
%           contains the triangle vertices.  These are the point indices of
%           the the three maxima that make up the triangle (the index is
%           into DT.X)
%
%        maxclass(i).geo_tris{:, 2}
%           contains the triangle vertices.  These are the maxclass indices of
%           the the three maxima that make up the triangle (the index is
%           into maxclass)
%
%        maxclass(i).geo_tris{:, 3}
%           contains a 1 by k vector which contains the sequence of
%           vertices making up the complete geodesic triangle.  This
%           includes all of the poins in each of the three geodesics as
%           well as the enpoints of each geodesic (the maxima).  This is
%           primarily for 2D plotting.  This allows us to color in the
%           geodesic triangles using a patch.  In 3D this won't work.
%
%        maxclass(i).geo_tris{:, 4}
%           contains 1 or 0 depending if the triangle is included or not
%           (respectively). A triangele is included if all three geodesics
%           making up the "sides" of the triangele were included by the
%           SelectGeodesics function.  In this case we put a 1 in the
%           fourth column.  Otherwise we put a 0 to indicate the triangle
%           is not included.
%        

Geodesic_Tris = cell(length(maxclass),3);

% First find all *possible* geodesic triangles.  This is equivalent to
% finding all triangles in a graph (three pairwise adjacent vertices). The
% graph here consists of the maxima as vertices and the geodesics as edges.
for m = 1:length(maxclass)
   NbsID = maxclass(m).nbormaxid;
   if isempty(NbsID)
      % then max m has no max neighbors so it can't be in a triangle.
      continue
   end
   
   % Transpose if needed for iterating the next for loop
   if size(NbsID,2) == 1 && size(NbsID,1) > 1
      NbsID = NbsID';
   end
   
   for a = NbsID
      NbsNbsID = maxclass(a).nbormaxid;
      CommonNbs = intersect(NbsID',NbsNbsID);
      if isempty(CommonNbs)
         continue
      end
      % Transpose if needed for iterating the next for loop
      if size(CommonNbs, 2) == 1 && size(CommonNbs, 1) > 1
         CommonNbs = CommonNbs';
      end
      for b = CommonNbs
         TriangleMaxID = sort([m, a, b]);
         TrianglePointID = sort([maxclass(m).max,....
            maxclass(a).max, maxclass(b).max]);
         Geodesic_Tris{m,1} = [Geodesic_Tris{m,1};TriangleMaxID];
         Geodesic_Tris{a,1} = [Geodesic_Tris{a,1};TriangleMaxID];
         Geodesic_Tris{b,1} = [Geodesic_Tris{b,1};TriangleMaxID];
         Geodesic_Tris{m,2} = [Geodesic_Tris{m,2};TrianglePointID];
         Geodesic_Tris{a,2} = [Geodesic_Tris{a,2};TrianglePointID];
         Geodesic_Tris{b,2} = [Geodesic_Tris{b,2};TrianglePointID];
      end
   end
end

% Remove redudant triangles.  Every triangle was found 6 times above.
for a = 1:length(Geodesic_Tris)
   [Geodesic_Tris{a,1}, r, ~] = unique(Geodesic_Tris{a},'rows');
   [Geodesic_Tris{a,2}] = Geodesic_Tris{a,2}(r,:);
end

% Allocate the cells in the 3rd column of Geodesic_Tris.  Also preallocate
% for maxclass(a).geo_tris (we will catalog the triangles in two places).
% Note that: length(Geodesic_Tris) == length(maxclass)
for a = 1:size(Geodesic_Tris, 1) % each max
   Geodesic_Tris{a,3}   = cell(size(Geodesic_Tris{a,1},1),2);
   maxclass(a).geo_tris = cell(size(Geodesic_Tris{a,1},1),4);
end

% Now the main block of the function. In this block we create the sequence
% of points which outline each geodesic triangle.  This is saved in the
% third column of Geodesic_Tris. This sequence is a list of indices to
% points in the raw data matrix, DT.X. This information is used for
% plotting later (for 2D only).  We also store each geodesic triangle in
% the maxclass struct. Also, we mark each triangle as being included or
% not.
for a = 1:length(Geodesic_Tris) % each max
   for b = 1:size(Geodesic_Tris{a,1},1) % each geodesic triangle
      
      % Get the current triangle expressed with both indices
      S = Geodesic_Tris{a,1}(b,:); % maxclass indices
      T = Geodesic_Tris{a,2}(b,:); % DT.X (data) indices
      
      % The following block greatly speeds up this nested for loop so we
      % are not duplicating work.
      if S(1) < a
         % Then we have already cataloged this triangle.  Move on to the
         % next triangle (the next iteration of the for loop indexed by b).
         continue
      end
      
      % Note that the index b is the row of the current triangle in the
      % current max.  But we need a row index for where this triangle is
      % for the other two vertices in the current triangle.  This will
      % allow us to store everything for all three vertices during this
      % iteration.  Then we don't have to repeat work later in the loop.
      %
      % Get the logical row index for the other two vertices
      r2 = ismember(Geodesic_Tris{S(2),1}, S, 'rows');
      r3 = ismember(Geodesic_Tris{S(3),1}, S, 'rows');
      
      % Store the current triangle T and S in maxclass for all three
      % vertices of the triangle S (S indexes into maxclass, T indexes into
      % the raw data matrix DT.X)
      maxclass(a).geo_tris{b,1} = T;
      maxclass(S(2)).geo_tris{r2,1} = T;
      maxclass(S(3)).geo_tris{r3,1} = T;
      maxclass(a).geo_tris{b,2} = S;
      maxclass(S(2)).geo_tris{r2,2} = S;
      maxclass(S(3)).geo_tris{r3,2} = S;
      
      % Now construct the sequence of vertices which make up the current
      % triangle. We will use this for 2D plotting.
      temp = vertcat(maxclass(S(1)).geodesics{:,2});
      r = temp(:,2) == T(2);
      Geo_temp1 = maxclass(S(1)).geodesics{r,1};
      Included_Temp1 = maxclass(S(1)).geodesics{r,4};
      
      temp = vertcat(maxclass(S(2)).geodesics{:,2});
      r = temp(:,2) == T(3);
      Geo_temp2 = maxclass(S(2)).geodesics{r,1};
      Included_Temp2 = maxclass(S(2)).geodesics{r,4};
      
      temp = vertcat(maxclass(S(3)).geodesics{:,2});
      r = temp(:,2) == T(1);
      Geo_temp3 = maxclass(S(3)).geodesics{r,1};
      Included_Temp3 = maxclass(S(3)).geodesics{r,4};
      
      % Store the geodesic triangle in the current triangle slot
      Geodesic_Tris{a,3}{b,1} = [Geo_temp1,...
         Geo_temp2(2:end),...
         Geo_temp3(2:end-1)];
      % Store this same triangle in the slot for the other 2 vertices.  We
      % use the logical indices found earlier.  Also store in maxclass.   
      Geodesic_Tris{S(2), 3}{r2, 1} = Geodesic_Tris{a,3}{b,1};
      Geodesic_Tris{S(3), 3}{r3, 1} = Geodesic_Tris{a,3}{b,1};
      maxclass(a).geo_tris{b,3} = Geodesic_Tris{a,3}{b,1};
      maxclass(S(2)).geo_tris{r2,3} = Geodesic_Tris{a,3}{b,1};
      maxclass(S(3)).geo_tris{r3,3} = Geodesic_Tris{a,3}{b,1};
      
      
      % Now mark whether every geodesic surrounding the current triangle
      % was selected.  If all geodesics were selected, then put a 1.  If
      % they weren't, then the triangle is not complete so put a 0.
      if Included_Temp1 == 2 && Included_Temp2 == 2 && Included_Temp3 == 2
         % Then the triangle is included.  Put a 1
         Geodesic_Tris{a,3}{b,2} = 1;
         Geodesic_Tris{S(2),3}{r2,2} = 1;
         Geodesic_Tris{S(3),3}{r3,2} = 1;
         maxclass(a).geo_tris{b,4} = 1;
         maxclass(S(2)).geo_tris{r2,4} = 1;
         maxclass(S(3)).geo_tris{r3,4} = 1;
      else
         % The triangle is not included.
         Geodesic_Tris{a,3}{b,2} = 0;
         Geodesic_Tris{S(2),3}{r2,2} = 0;
         Geodesic_Tris{S(3),3}{r3,2} = 0;
         maxclass(a).geo_tris{b,4} = 0;
         maxclass(S(2)).geo_tris{r2,4} = 0;
         maxclass(S(3)).geo_tris{r3,4} = 0;
      end
   end
end

end %function