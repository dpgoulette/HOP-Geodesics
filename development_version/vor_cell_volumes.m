function  Vols=vor_cell_volumes(VV,VC,GoodIndex)
% volumes - finds the volumes of the delaunay tetrahedra. Note that an
%       important part of this script relies on the fact that the point at
%       infinity is stored in the first entry in the cell.  
%
%       V is the matrix of voronoi vertices
%       C is the cell containing the indexes of each V-cel

% create a volume list as long as the number of data points.  Bad data
% points will ultimately have a 0 for its density.
Vols=zeros(length(VC),1);

for a=1:length(GoodIndex)
   % get the V-cell verts
   VcellVerts = VV( VC{GoodIndex(a)}, :);
   [~, Vols(GoodIndex(a))]=convhull(VcellVerts); %volume of V-cell
    
    % Update progress to the terminal.
    if mod(a,100)==0 && mod(a,1000) ~= 0
       fprintf('.')
    end
    if mod(a,1000)==0
        fprintf('\nCompleted %d Voronoi Cell Volumes\n',a);
    end
end
fprintf('\n')
end%function
