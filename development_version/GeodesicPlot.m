function GeodesicPlot(DT, GoodIndex, maxclass, maxindex)
% GeodesicPlot -- creates a plot of the HOP geodesics model.

fprintf('We will now plot the max geodesics.\n\n')
while true
   fprintf('Pause between maxima?\n')
   fprintf('  1) Yes\n')
   fprintf('  2) No\n')
   pause_option = input('(Enter 1 or 2): ');
   if pause_option == 1 || pause_option == 2
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end
fprintf('\n')

if size(DT.X,2) == 3
   GeosOnly3D(DT, GoodIndex, maxclass, maxindex, pause_option)
   
else % The data is 2 dimensional
   
   switch pause_option
      case 1 %pause between maxima
         
         fprintf('Would you like to see the step functions for each max?\n')
         while true
            fprintf('  1) Yes\n')
            fprintf('  2) No\n')
            plot_option = input('(Enter 1 or 2): ');
            if plot_option==1 || plot_option==2
               break
            else
               fprintf('\nERROR: You must enter 1 or 2.\n')
            end
         end
         fprintf('\n')
         
         if plot_option == 1
            GeosWithStepFunc(DT,GoodIndex, maxclass,maxindex)
         else
            GeosOnly2D(DT, GoodIndex, maxclass, maxindex, 1)
         end
         
      otherwise % No pause between maxima
         fprintf('\n')
         GeosOnly2D(DT, GoodIndex, maxclass, maxindex, pause_option)
   end
end

end % main function: GeodesicPlot

