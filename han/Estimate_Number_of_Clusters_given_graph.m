function [K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W, NUMC)
%%%This function estimates the number of clusters given the two huristics
%%%given in the supplementary materials of our nature method paper
%W is the similarity graph
%NUMC is a vector which contains the possible choices of number of
%clusters.


%%K1 is the estimated best number of clusters according to eigen-gaps
%%K12 is the estimated SECOND best number of clusters according to eigen-gaps

%%K2 is the estimated number of clusters according to rotation cost
%%K22 is the estimated SECOND number of clusters according to rotation cost

% an example would be [K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W,
% [2:5]);

%%Note that this function can only give an estimate of the number of
%%clusters. How to determine the "OPTIMAL" number of clusters, is still an
%%open question so far. 





if nargin < 2
    NUMC = 2:5;
end

if min(NUMC)==1
    warning('Note that we always assume there are more than one cluster.');
    NUMC(NUMC<=1) = [];
end
W = (W + W')/2;
W = W - diag(diag(W));

if ~isempty(NUMC)
    degs = sum(W, 2);
    D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
    % compute unnormalized Laplacian
    L = D - W;
    degs(degs == 0) = eps;
    % calculate D^(-1/2)
    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
    % calculate normalized Laplacian
    L = D * L * D;
    % compute the eigenvectors corresponding to the k smallest
    % eigenvalues
    [U, eigenvalue] = eigs(L, max(NUMC)+1, eps);
    eigenvalue  = diag(eigenvalue);
    [a,b] = sort((eigenvalue),'ascend');
    eigenvalue = (eigenvalue(b));
    U = U(:,b);
    eigengap = abs(diff(eigenvalue));
    eigengap = eigengap.*(1-eigenvalue(1:end-1))./(1-eigenvalue(2:end));
    for ck = NUMC
        Cindex = find(NUMC==ck);
        UU = U(:,1:ck);
        %UU = UU./repmat(sqrt(sum(UU.^2,2)),1,size(UU,2));
        [EigenvectorsDiscrete,EigenVectors ]=discretisation(UU);
        EigenVectors = EigenvectorsDiscrete.^2;
        [temp1,temp] = sort(EigenVectors,2, 'descend');
        quality(Cindex) = (1-eigenvalue(ck+1))/(1-eigenvalue(ck))*sum(sum(diag(1./(temp1(:,1)+eps))*temp1(:,1:max(2,ck-1))));
    end
    
    [tt1, t1] = sort(eigengap(NUMC),'descend');K1 = NUMC(t1(1));K12 = NUMC(t1(2));
    [tt2, t2] = sort(quality,'ascend');K2 = NUMC(t2(1));K22 = NUMC(t2(2));
   
end


