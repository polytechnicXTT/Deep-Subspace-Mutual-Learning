function [ Y ] = simpleLP( W,Y, T)
if nargin<3
    T  = 20;
end
W = dn(W,'ave');
Y0 = Y;
% 
labelindex = find(sum(Y==0,2)==(size(Y,2)-1));
testindex = setdiff([1:length(W)], labelindex);
Y0(testindex,:)=0;

% Yt = Y(testindex,:);
% 
% [a,b] = sort(Yt,2, 'descend');
% 
% index=find(abs(Yt(sub2ind(size(Yt),[1:size(Yt,1)]',b(:,1)))./(Yt(sub2ind(size(Yt),[1:size(Yt,1)]',b(:,2)))+eps))<1.2);
% Yt(index,:)=0;
% 
% index=find(abs(Yt(sub2ind(size(Yt),[1:size(Yt,1)]',b(:,1)))./(Yt(sub2ind(size(Yt),[1:size(Yt,1)]',b(:,2)))+eps))>3);
% if ~isempty(index)
%     Yt(index,:)=0;
%     Yt(sub2ind(size(Yt),[index],b(index,1))) = 1;
%     Y(testindex,:) = Yt;
% end
% 
% 
% labelindex = find(sum(Y==0,2)==(size(Y,2)-1));

%Y = (Y+ Y0)/2;
Y = Y0;
for i = 1:T
    Y = W*Y;
    Y(labelindex,:) = Y0(labelindex,:);
end
%
% Y = 0.01*(eye(length(W))-0.99*W)\Y;
%
% Y(labelindex,:) = Y0(labelindex,:);
end



