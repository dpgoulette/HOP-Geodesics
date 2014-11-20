%assume r is a reflexive relation, represented rowwise
function shape=constructshape(r,d)
	%recreate reflexive relation
    shape=[];
	e=sortrows(unique(sort(r,2),'rows'),1);
	e=sortrows([e;e(:,[2 1])],1);
	for i=reshape(unique(e(:,1)),1,length(unique(e(:,1))))
        s=e(e(:,1)==i,:);
        if (~isempty(s))
            for j=2:(d+1)
                t=[];
                for k=1:size(s,1)
                    v=e(e(:,1)==s(k,j),:);
                    t=[t;[repmat(s(k,1:(j-1)),size(v,1),1) v]];
                end
                s=t;
            end
            s=s(s(:,1)==s(:,end),1:(end-1));
            %s=sortrows(unique(sort(s,2),'rows'),1);
            s=s(~sum(diff(sort(s,2)',1)'==0,2),:);
            %eliminate the circulant/reverse circulant ones
            j=1;
            while (j<=size(s,1))
                p=gallery('circul',s(j,:));
                p=[p;fliplr(p)];
                s=s(boolean([repmat(1,j,1);~ismember(s((j+1):end,:),p,'rows')]),:);
                j=j+1;
            end
            shape=[shape;s(~sum(diff(sort(s,2)',1)'==0,2),:)];
        end
        e=e(~sum(e==i,2),:);
	end
end