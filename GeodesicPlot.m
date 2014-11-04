
fprintf('This script plots the results of HOP and max geodesics.\n')
fprintf('Choose one of the following options: ')
while true
   fprintf('  1) Plot the raw data, geodesics AND the geodesic step functions.\n');
   fprintf('  2) Plot ONLY the raw data and geodesics.\n')
   plot_option = input('Which of the above options do you prefer: ');
   if plot_option==1 || plot_option==2 
      break
   else
      fprintf('\nERROR: You must enter 1 or 2.\n')
   end
end

if plot_option == 1
   clear plot_option
   GeosWithStepFunc(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
else
   clear plot_option
   GeosOnly(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
end
