
fprintf('This script plots the results of HOP and max geodesics.\n')
fprintf('Would you like to see the geodesics only? Or do you want \n')
fprintf('to see the step functions for each max as well?\n')
while true
   fprintf('  1) Geodesics only.\n')
   fprintf('  2) Geodesics with step functions.\n')
   plot_option = input('Which of the above options do you prefer: ');
   if plot_option==1 || plot_option==2 
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end
fprintf('\n')
fprintf('\n')
pause(.2)
if plot_option == 1
   clear plot_option
   GeosOnly(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
else
   clear plot_option
   GeosWithStepFunc(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
end
