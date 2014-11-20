%assume r is a reflexive relation, represented rowwise
function simplex=constructsimplex(r,d)
	%recreate reflexive relation
    simplex=[];
	e=sortrows(unique(sort(r,2),'rows'),1);
	e=sortrows([e;e(:,[2 1])],1);
	for i=reshape(unique(e(:,1)),1,length(unique(e(:,1))))
        n=e(e(:,1)==i,2);
        if (length(n)>=d)
            ne=sort(nchoosek(n,d),2);
            nem=zeros(size(ne,1),1);
            for j=1:size(ne,1)
                if (nchoosek(size(ne,2),2)==sum(ismember(sort(nchoosek(ne(j,:),2),2),e,'rows')))
                    nem(j)=1;
                end
            end
            ne=ne(boolean(nem),:);
            simplex=[simplex;[repmat(i,size(ne,1),1),ne]];
            e=e(~sum(e==i,2),:);
        end
	end
end
