function GoodEdges = AlphaCellsSelect2D(DT,GoodEdges,VV,VC)
%

% Already have good delaunay edges (1-cells).  Calculate the value of
% alpha for each good edge
fprintf('Calculating alpha for the 1-cells.\n\n')
cells1 = EpsilonOneCells2d(DT,GoodEdges,VV,VC);

%       %%%%  THIS BLOCK IS DISABLED %%%%%%
%       % If we want to use the 2-cells at some time enable this block and
%       % make sure that the function returns the cells2 matrix
%
%       % We already have the 1-cells (GoodEdges).  We need the 2-cells
%       DTtris=DT.Triangulation; % 2-cells
%       % Now remove bad triangles.  They are "bad" if they contain
%
%       [R,~] = find(ismember(DTtris, BadDataID));
%       R = unique(R);
%       cells2 = DTtris;
%       cells2(R,:) = [];
%
%       % Find the circumcenter radius for the triangles.  This is the value of
%       % alpha for all of teh 2 cells.
%       [~,RCC] = circumcenters(DT);
%       % Now delete the radii from bad triangles
%       RCC(R)=[];
%       % Save the 2 cells along with their alpha value in the third column.
%       cells2=[cells2, RCC];
%
%%%%% WE WOULD ALSO NEED SOME SORT OF A SELECTION SCHEME!!! RIGHT??!!

% User selects which alpha complex selection scheme they want.
fprintf('Which alpha 1-cell selection scheme do you want?\n')
while true
   fprintf('   1) Choose a global alpha for all edges (as a percentage).\n')
   fprintf('   2) Standard deviations above the mean edge alpha.\n')
   fprintf('   3) Our parameter free local alpha select.\n')
   alpha_select_option = input('Choose one of the above: ');
   if alpha_select_option == 2 || alpha_select_option == 1 ||...
         alpha_select_option == 3
      break
   else
      fprintf('ERROR! You must enter 1 or 2.\n\n')
   end
end
fprintf('\n')
switch alpha_select_option
   case 1
      while true
         fprintf('\nWhat percent of the smallest edge alpha values do you want ')
         fprintf('to keep? \n')
         keep_percent = input('Enter a decimal from 0 to 1: ');
         if keep_percent < 0 || keep_percent > 1
            fprintf('\nERROR: The number must be in the interval [0,1].\n')
         else
            break
         end
      end
      cells1 = sortrows(cells1,3);
      num_keep = floor(size(cells1,1) * keep_percent);
      if num_keep <= 0
         fprintf('\nAll edges will be removed!!\n')
         GoodEdges =[];
      elseif num_keep >= size(cells1,1)
         fprintf('\nNo edges will be deleted. This is the full Delaunay.\n')
      else
         cells1(num_keep + 1:end,:) = [];
         GoodEdges = cells1(:,[1,2]);
      end
      
   case 2
      fprintf('\nHow many standard deviations above the mean do you want')
      fprintf(' to keep? Entering a negative number will be below the mean.')
      fprintf(' Every edge with alpha less than or equal to the number of')  
      fprintf(' standard deviations will be kept.\n')
      standard_devs = input('   Enter any real number: ');
      
      cells1 = sortrows(cells1,3);
      SD = std(cells1(:,3));
      cutoff = mean(cells1(:,3)) + standard_devs * SD;
      cutoff_start = find(cells1(:,3) > cutoff, 1);
      if cutoff_start <= 0
         fprintf('\nAll edges will be removed!!\n')
         GoodEdges =[];
      elseif cutoff_start >= size(cells1,1)
         fprintf('\nNo edges will be deleted. This is the full Delaunay.\n')
      else
         cells1(cutoff_start:end,:) = [];
         GoodEdges = cells1(:,[1,2]);
      end
   otherwise
      error('\n\nNot yet implemented\n\n')
end
end