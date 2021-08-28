function [W]=SNF(Wall,K,t,ALPHA)

if nargin < 2
    K = 20;
end
if nargin < 3
    t = 20;
end

if nargin < 4
    ALPHA = 1;
end

C = length(Wall);
[m,n]=size(Wall{1});
for i = 1 : C
    Wall{i} = Wall{i}./repmat(sum(Wall{i},2),1,n);
    Wall{i} = (Wall{i} + Wall{i}')/2;
end

for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end

Wsum = zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
for ITER=1:t
    for i = 1 : C
        %Wall0{i}=newW{i}*(0.95*(Wsum - Wall{i})/(C-1)+0.05*eye(length(Wsum)))*newW{i}';
        Wall0{i}=newW{i}*(Wsum - Wall{i})*newW{i}'/(C-1);
    end
    for i = 1 : C
        Wall{i} = BOnormalized(Wall0{i},ALPHA);
    end
    Wsum = zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
    end
%     
end

W = Wsum/C;
W = W./repmat(sum(W,2),1,n);;
W = (W +W'+eye(n))/2;
end

function W = BOnormalized(W,ALPHA)
if nargin < 2
    ALPHA = 1;
end
W = W+ALPHA*eye(length(W));
W = (W +W')/2;
end

function newW = FindDominateSet(W,K)
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
clear IW2;
clear temp;
end

