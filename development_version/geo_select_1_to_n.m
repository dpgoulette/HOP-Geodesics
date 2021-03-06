function maxclass = geo_select_1_to_n(maxclass, total_num_geodesics)
%  geo_select_1_to_n -- loops through each max in maxclass and determines which
%        geodesics to keep and which to throw away.  We rank all of the
%        geodesics in order of their length (based on our metric).  The
%        rank of the shortest geodesic in ALL of maxclass has rank 1 and
%        the longest geodesic has rank r = "the total number of geodesics."
%        The selection scheme for the i-th max is based on the rank order
%        of the lengths of the geodesics connected to max i.  We use a
%        persistence-like scheme to find out the cut-off point for which
%        geodesics are connected to the i-th max and which are removed. For
%        example, if all of the geodesics connected to i are very short (in
%        comparison to all of the geodesics in maxclass), then this max is
%        close to it's neighbor max. So it is in a denser region, and hence
%        we keep all of it's geodesics.  But suppose, for a different
%        example, max i has 6 geodesics connected to it.  Further suppose
%        that the two shortest geodesics are very short, but the longest 4
%        are very long, then there will be a large gap in the rank between
%        the 2nd geodesic and the 3rd (when taken in order of length).  Our
%        algorithm finds the largest gap in the ranks and keeps all of the
%        geodesics before that gap.  So in this example we would keep the 2
%        shortest and not include the 4 longest.  NOTE! Since there are n
%        geodesics connected to max i, there are n+1 gaps to consider. This
%        is because we compare the rank of the shortest geodesic to 1 and
%        we compare the rank of the longest geodesic to r = "the total
%        number of geodesics.  HOWEVER, in this function we do not consider
%        gap zero (the difference in rank between the shortest geodesic
%        connected to i and the shortest geodesic in all of maxclass).
%        Thus, this function will include AT LEAST ONE GEODESIC FOR EVERY
%        MAX IN MAXCLASS.  After this function terminates, the parent
%        function, SelectGeodesics, will only include geodesics that were
%        selected from BOTH ends (meaning the maxima at both ends of the
%        geodesic selected the geodesic connecting them).
%
%        geo_select_1_to_n updates the fourth column in maxclass(i).geodescis for
%        each geodesic connected to i.


for Max = 1:length(maxclass)
   % There are cases where a maxclass has no geodesics connected to
   % it.  This only happens on the boundary of the data set or if you
   % hop on an alpha complex. Skip to the next max if this is the
   % case.
   if isempty(maxclass(Max).nbormaxid)
      continue
   end
   NumMaxNeighbors = length(maxclass(Max).nbormaxid);
   % This block selects the longest of bars 1 through n (not 0).
   % So this max must select at least the shortest geodesic
   % attached to it.
   temp = zeros(NumMaxNeighbors, 1);
   % Note that the for loop starts at 1 since we are not considering
   % the zero bar. (Note: this means that every max will select at
   % least one geodesic.  This is o.k. because, in the end, we will
   % only include geodesics that are selected from BOTH ends.)
   for Bar = 1 : NumMaxNeighbors
      if Bar < NumMaxNeighbors %not the first bar AND not the last
         %bar length is the rank of the next geodesic minus the
         %rank of the current geodesic.
         temp(Bar) = maxclass(Max).geodesics{Bar + 1, 3} - ...
            maxclass(Max).geodesics{Bar, 3};
      else %the last bar.  So Nbor == NumMaxNeighbors +1
         %bar length is the rank of the longest geodesic minus the
         %rank of the current geodesic.
         temp(Bar) = total_num_geodesics - ...
            maxclass(Max).geodesics{Bar, 3};
      end
   end
   [~,slot] = max(temp);
   NumGeosSelected = slot;
   
   % Now put a 1 in the fourth column for each geodesic that was
   % selected in the above process.  Put 0 for the rest that weren't
   % selected.
   iter = 1;
   while iter <= NumGeosSelected
      maxclass(Max).geodesics{iter, 4} = 1;
      iter = iter + 1;
   end
   while iter <= NumMaxNeighbors
      maxclass(Max).geodesics{iter, 4} = 0;
      iter = iter + 1;
   end
   
end %loop through maxclass

end % function: geo_select_1_to_n