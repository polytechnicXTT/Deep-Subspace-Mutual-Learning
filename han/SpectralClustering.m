function [group, eigengap] = SpectralClustering(W, NUMC)
%SPECTRALCLUSTERING Executes spectral clustering algorithm

% calculate degree matrix
degs = sum(W, 2);
D  = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
L = D - W;
k = max(NUMC);
% compute normalized Laplacian if needed

% avoid dividing by zero
degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
% calculate normalized Laplacian
L = D * L * D;

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
[U, eigenvalue] = eigs(L, k, eps);
[a,b] = sort(diag(eigenvalue),'ascend');
eigenvalue = eigenvalue(:,b);
U = U(:,b);
eigengap = abs(diff(diag(eigenvalue)));
U = U(:,1:k);
% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
%U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
%U = U./repmat(sqrt(sum(U.^2,2)),1,size(U,2));

flag =0;
for ck = NUMC
    Cindex = find(NUMC==ck);
    UU = U(:,1:ck);
    UU = UU./repmat(sqrt(sum(UU.^2,2)),1,size(UU,2));
    [EigenvectorsDiscrete]=discretisation(UU);
    [~,temp] = max(EigenvectorsDiscrete,[],2);
%     for i = 1 : ck
%         initcenter(i,:) = mean(UU(temp==i,:));
%     end
    
    Cluster{Cindex} = temp;
end


if length(NUMC)==1
    group=Cluster{1};
else
    group = Cluster;
end


end