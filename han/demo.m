clc
clear
close all
%%%Load the data

%load simulation.mat

%Data1 is of size n x d_1, where n is the number of patients, d_1 is the number of genes, e.g.
%Data2 is of size n x d_2, where n is the number of patients, d_2 is the number of methylation, e.g.

%%%First, set all the parameters.
K = 20;%number of neighbors, usually (10~30)
alpha = 0.4; %hyperparameter, usually (0.3~0.8)0.5
T = 20; %Number of Iterations, usually (10~20)

%If the data are all continuous values, we recommend the users to perform standard normalization before using SNF, though it is optional depending on the data the users want to use. 

%Data1 = data1;Data2 = data2;
Data1 = xlsread('E:\PyCharm Project\train\maybe_best\W7_0.xlsx','B2:JR278')';
Data2 = xlsread('E:\PyCharm Project\train\maybe_best\W7_1.xlsx','B2:JR278')';
Data3 = xlsread('E:\PyCharm Project\train\maybe_best\W7_2.xlsx','B2:JR278')';

Data1 = Standard_Normalization(Data1);
Data2 = Standard_Normalization(Data2);
Data3 = Standard_Normalization(Data3);

%%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
%if the data is discrete, we recommend the users to use chi-square distance
Dist1 = dist2(Data1,Data1);
Dist2 = dist2(Data2,Data2);
Dist3 = dist2(Data3,Data3);

%%%next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha);
W2 = affinityMatrix(Dist2, K, alpha);
W3 = affinityMatrix(Dist3, K, alpha);

%%% These similarity graphs have complementary information about clusters.
%displayClusters(W1,label);
%displayClusters(W2,label);

%next, we fuse all the graphs
% then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF({W1,W2,W3}, K, T);
%%%%With this unified graph W of size n x n, you can do either spectral clustering or Kernel NMF. If you need help with further clustering, please let us know. 
%%for example, spectral clustering
C = 3;%%%number of clusters5
group = SpectralClustering(W,C);%%%the final subtypes information

%%%you can evaluate the goodness of the obtained clustering results by calculate Normalized mutual information (NMI): if NMI is close to 1, it indicates that the obtained clustering is very close to the "true" cluster information; if NMI is close to 0, it indicates the obtained clustering is not similar to the "true" cluster information.

displayClusters(W,group);

SNFNMI = Cal_NMI(group, label);


%%%you can also find the concordance between each individual network and the fused network
ConcordanceMatrix = Concordance_Network_NMI({W,W1,W2},C);


%%%Here we provide two ways to estimate the number of clusters. Note that,
%%%these two methods cannot guarantee the accuracy of esstimated number of
%%%clusters, but just to offer two insights about the datasets.

[K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W, [2:5]);

fprintf('The best number of clusters according to eigengap is %d\n', K1);
fprintf('The best number of clusters according to rotation cost is %d\n', K2);



%%%%We also provide an example using label propagation to predict the labels of new data points(see Demo2.m).


