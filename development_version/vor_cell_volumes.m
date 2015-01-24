function  Vols=vcellvolumes(V,C)
% volumes - finds the volumes of the delaunay tetrahedra. Note that an
%       important part of this script relies on the fact that the point at
%       infinity is stored in the first entry in the cell.  
%
%       V is the matrix of voronoi vertices
%       C is the cell containing the indexes of each V-cel

Vols=zeros(length(C),1);

for a=1:length(C)
    if C{a}~=1 %check that the cell does NOT contain point at infinity
        VcellVerts=V(C{a},:);%use the C{a} index vector to get V-cell verts
        [~, Vols(a)]=convhull(VcellVerts); %volume of V-cell
    end
    
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
