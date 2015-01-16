function maxclass = geodesic_classify(maxclass)
% adds a 5th column to the cell in maxclass(i).geodesics.  This column will
% hold a 0, 1, 2 or 3.
%     0 if the geodesic is not included
%     1 if the geodesic is included but not in a triangle
%     2 if the geodesic is included, in a triangle but not in a tet
%     3 if teh geodesic is included in a tet.
%
% NOTE!!  This function only includes a tetra if it was selected by both
% maxima at its endpoints.  That is, only if there is a 2 in
% maxclass(i).geodesics{j,4}.   If there is a 1 in that slot, then the
% geodesic was only selected from one end.  So it is not included.
% 

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
   
   % Now check if there are any tetras.     
   if ~isempty(maxclass(a).geo_tetras) &&...
         any(vertcat(maxclass(a).geo_tetras{:,3}));
      % Then there are tetras.  And at least one of them is included.
      % Find all of the geodesics which are included in tetras, and include
      % the current max point.  They need to be expressed in the order that
      % they are in maxclass(i).geodesics{j,2}(:,[1,2])
      geos = find_tet_geos(geos, maxclass, a);
      if any(isnan(geos(:,3)))
         geos = find_tri_geos(geos, maxclass, a);
      end
   elseif ~isempty(maxclass(a).geo_tris) &&...
         any(vertcat(maxclass(a).geo_tris{:,4}));
      % Thus there were no tetras included but there are triangles included
      geos = find_tri_geos(geos, maxclass, a);
   end
   
   if any(isnan(geos(:,3)))
      geos(isnan(geos)) = 1;
   end
   
   maxclass(a).geodesics(:,5) = num2cell(geos(:,3));
   
end
end % main function

function geos = find_tet_geos(geos, maxclass, a)
%

% get the current max's tetras
temp = vertcat(maxclass(a).geo_tetras{:,1});
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
%

% get the current max's tris
temp = vertcat(maxclass(a).geo_tris{:,1});
% get the unique vertices
temp = unique(temp);
% remove the current max.  we want the neighbors
temp(temp == maxclass(a).max) = [];
% find the geos that are in tetras
r = ismember(geos(:,2), temp) & isnan(geos(:,3));
% assign 3 to those rows.  They are 3D geos
geos(r,3) = 2;
end % find_tri_geos