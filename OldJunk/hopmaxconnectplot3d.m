%voronoi(DT,'g')
plot3(DT.X(GoodIndex,1),DT.X(GoodIndex,2),DT.X(GoodIndex,3),'g.','markersize',5)

%plot3(DT.X(BadDataID,1),DT.X(BadDataID,2),DT.X(BadDataID,3),'y.','markersize',1)
hold on;
axis equal;

Mid=vertcat(maxclass.max);
plot3(DT.X(Mid,1),DT.X(Mid,2),DT.X(Mid,3),'b*','markersize',5)
%plot3(DT.X(Mid,1),DT.X(Mid,2),DT.X(Mid,3),'k.','markersize',5)

threshold = 1;
pause

%plot(DT.X(GoodIndex,1),DT.X(GoodIndex,2),'g.')

%plot(DT.X(BadDataID,1),DT.X(BadDataID,2),'r.')
% Mid=vertcat(maxclass.max);
% plot3(DT.X(Mid,1),DT.X(Mid,2),DT.X(Mid,3),'b.','markersize',3)
%plot(DT.X(Mid,1),DT.X(Mid,2),'k+')
A=maxconnect;
z=vertcat(A{:,1});
z(:,[1 2])=[];
[~,id]=sort(z);
maxconnectsorted=A(id,:);

for a=1:size(maxconnectsorted,1)
    if maxconnectsorted{a,1}(3) < threshold
        E=zeros(length(maxconnectsorted{a,2})-1,2);
        for b=1:size(E,1)
            E(b,:)=[maxconnectsorted{a,2}(b),maxconnectsorted{a,2}(b+1)];
        end
        X=[DT.X(E(:,1),1)';DT.X(E(:,2),1)'];
        Y=[DT.X(E(:,1),2)';DT.X(E(:,2),2)'];
        Z=[DT.X(E(:,1),3)';DT.X(E(:,2),3)'];
        
%         [~,Xb] = find((X < -.08 | X > -0.04));
%         X(:,Xb)=[];
%         Y(:,Xb)=[];
%         Z(:,Xb)=[];
%         [~,Yb] = find(Y < -.02);
%         X(:,Yb) =[];
%         Y(:,Yb) =[];
%         Z(:,Yb) =[];

        plot3(X,Y,Z,'r-')
        if mod(a,50)==0
            title(['Edge length threshold ',...
                num2str(maxconnectsorted{a,1}(3))])
            pause(0.2)
        end
        if mod(a,50)==0
            pause
        end
    else
        break
    end
end
clear a b z Mid X Y z id A
clear E
