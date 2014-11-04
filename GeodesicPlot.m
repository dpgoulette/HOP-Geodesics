
fprintf('We will now plot the max geodesics.\n\n')
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
fprintf('\n')
pause(.2)
if plot_option == 1
   clear plot_option
   GeosWithStepFunc(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
else
   clear plot_option
   GeosOnly(DT,GoodIndex,maxconnect,maxclass,maxindex,GoodMaxGeodesics)
end
