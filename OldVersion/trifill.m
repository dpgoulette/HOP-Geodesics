% figure
% P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.');
% axis equal
% hold on
% axis tight
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% P2 = plot(DT.X(maxindex,1),DT.X(maxindex,2),'g.'); 


% z=[vertcat(maxconnect{:,1}),vertcat(maxconnect{:,4})];
% z(:,3)=[];
% z(z(:,3)~=2,:)=[];
% z(:,3)=[];
% 
% X=[DT.X(z(:,1),1)';DT.X(z(:,2),1)'];
% Y=[DT.X(z(:,1),2)';DT.X(z(:,2),2)'];
% E = plot(X,Y,'m-');
% 
% geo_tris = constructsimplex(z,2);
% 
% tri_list = cell(size(geo_tris,1),1);
% for i = 1:size(tri_list,1)
%    tri_list{i} = [geo_tris(i,1), geo_tris(i,2);...
%                   geo_tris(i,2), geo_tris(i,3);...
%                   geo_tris(i,3), geo_tris(i,1);];
% end
% 
% edge_plot_list = cell(length(tri_list),1);
% fill_plot_list = cell(length(tri_list),1);
% for i=1:size(tri_list,1)
%    X=[DT.X(tri_list{i}(:,1),1)';DT.X(tri_list{i}(:,2),1)'];
%    Y=[DT.X(tri_list{i}(:,1),2)';DT.X(tri_list{i}(:,2),2)'];
% %    edge_plot_list{i} = plot(X,Y,'m-');
%    x = X(1,:);
%    y = Y(1,:);
%    fill_plot_list{i} = patch(x,y,'magenta');
% end
PlotHandles = cell(size(Geodesic_Tris,1),1);
for a = 1:size(Geodesic_Tris,1)
   for b = 1:size(Geodesic_Tris{a,3},1)
      if a == Geodesic_Tris{a,1}(b,1)
         %Then we haven't plotted it yet.
      if Geodesic_Tris{a,3}{b,2} == 1
         geo_tri_temp = Geodesic_Tris{a,3}{b,1};
         x=DT.X(geo_tri_temp,1)';
         y=DT.X(geo_tri_temp,2)';
         PlotHandles{a,b} = patch(x,y,'magenta','FaceAlpha',.7,...
                                  'EdgeColor','magenta');
      end
      end
   end
end