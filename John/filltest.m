clf
P1 = plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'k.');
axis equal
hold on
axis tight
set(gca, 'XTick', []);
set(gca, 'YTick', []);
P2 = plot(DT.X(maxindex,1),DT.X(maxindex,2),'g.');
axis( [-0.000205496214662   0.005191413740226   0.097779668761888   0.102814262824578]);
P = cell(3,1);
x=[];
y=[];
for j = 1:3
   GeoPath = z{j,2};
   GeoEdges = [GeoPath(1:end-1)', GeoPath(2:end)'];
   X=[DT.X(GeoEdges(:,1),1)';DT.X(GeoEdges(:,2),1)'];
   Y=[DT.X(GeoEdges(:,1),2)';DT.X(GeoEdges(:,2),2)'];
   x=[x,X(1:end)];
   y=[y,Y(1:end)];
   P{j} = plot(X,Y,'r-');
end
patch(x,y,'y')


