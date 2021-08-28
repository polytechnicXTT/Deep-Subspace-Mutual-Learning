function [ave_s] = New_silhouette(X, clust)


D=affinityMatrix(X,20,0.5);

r=size(D,1);
s=zeros(r,1);
a=zeros(r,1);
b=zeros(r,1);
for i=1:1:r
    
ID1=find(clust==clust(i));
% a(i)=mean(D(i,ID1));
if (length(ID1)==1)
    a(i)=0;
else
a(i)=mean(D(i,ID1));
ID1(find(ID1==i))=[];
end

ID2=find(clust~=clust(i));
b(i)=max(D(i,ID2));
if max(a(i),b(i))==0
    s(i)=0;
else
s(i)=(b(i)-a(i))/(max(a(i),b(i)));
end
end
ave_s=mean(s);
end

%-------------------------------------------------------------------
