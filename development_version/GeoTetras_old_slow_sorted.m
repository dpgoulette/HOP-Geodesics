function tetras = GeoTetras_old_slow_sorted(maxclass)
%
tetras = cell(length(maxclass),3);
for a = 1:length(maxclass)
   if ~isempty(maxclass(a).geo_tris)
      for b = 1:size(maxclass(a).geo_tris, 1)
         current_tri = maxclass(a).geo_tris{b, 2};
         % Find any neighbors the 3 vertices have in common
         nbs_1 = maxclass(current_tri(1)).nbormaxid;
         nbs_2 = maxclass(current_tri(2)).nbormaxid;
         nbs_3 = maxclass(current_tri(3)).nbormaxid;
         common_nbs = sort(intersect(intersect(nbs_1, nbs_2), nbs_3));
         
         % Each common neighbor defines a tetra.  Store them all in the
         % Tetras cell for now.
         if ~isempty(common_nbs)
            for c = 1:length(common_nbs)
               tet_maxclass_id = sort([current_tri, common_nbs(c)]);
               tet_maxpoint_id = horzcat(maxclass(tet_maxclass_id).max);
               tetras{a,1} = [tetras{a,1}; tet_maxpoint_id];
               tetras{a,2} = [tetras{a,2}; tet_maxclass_id];
            end
         end
      end
   end
end

% Remove the duplicate Tetras
for a = 1:length(tetras)
   if ~isempty(tetras{a,1})
      tetras{a,1} = unique(tetras{a,1},'rows');
      tetras{a,2} = unique(tetras{a,2},'rows');
   end
end



end % main function

