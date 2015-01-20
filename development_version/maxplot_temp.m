function maxplot_temp(Geodesic_Tris,maxclass,DT,Max)

fprintf('\nPlotting max %d.\n',Max)

P1 = plot(DT.X(maxclass(Max).max,1),DT.X(maxclass(Max).max,2),'k*');
if ~isempty(Geodesic_Tris{Max})
   T = Geodesic_Tris{Max};
   T = unique(T);
   T(T==maxclass(Max).max)=[];
   Plot = cell(size(T));
   for a = 1:length(T)
      Plot{a} = plot(DT.X(T(a),1),DT.X(T(a),2),'m*');
   end
   pause
   
   delete(Plot{:})
end
delete(P1)
end