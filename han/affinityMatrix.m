function [W]=affinityMatrix(Diff,K,sigma);
% computes an affinity matrix for a given distance matrix
if nargin<3
    sigma = 0.5;
end
if nargin<2
    K = 20;
end

Diff=(Diff+Diff')/2;
Diff = Diff - diag(diag(Diff));
%Diff = sqrt(Diff+eps);
[T,INDEX]=sort(Diff,2);
[m,n]=size(Diff);
W=zeros(m,n);
TT=mean(T(:,2:K+1),2)+eps;
Sig=(repmat(TT,1,n)+repmat(TT',n,1) + 1*Diff)/3;
Sig=Sig.*(Sig>eps)+eps;
W=normpdf(Diff,0,sigma*Sig);

W = (W + W')/2;

return
