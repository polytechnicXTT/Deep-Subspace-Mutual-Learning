function ConcordanceMatrix = Concordance_Network_NMI(Wall,C)
%%This function outputs the concordance between all the networks.
%%SO if you have n networks, then Wall is a cell with n+1 networks among
%%which the first is the fused network
%%C is the number of cluster you want to cluster

%%the output ConcordanceMatrix is a (n+1)x(n+1) square matrix
N = length(Wall);
for i=1:N
    group{i} = SpectralClustering(Wall{i},C);
end
ConcordanceMatrix = zeros(N);
for i = 1:N-1
    for j = i+1:N
        ConcordanceMatrix(i,j) = Cal_NMI(group{i},group{j});
    end
end
ConcordanceMatrix = ConcordanceMatrix+ConcordanceMatrix'+eye(N);



end