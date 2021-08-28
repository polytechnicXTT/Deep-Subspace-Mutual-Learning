function Predicted_labels = group_predict(train,test,trainlabel,K, alpha, T)
%This function implements label propagation with multiple data types to
%predict the labels for new data.
%train is a mx1 cell in each of which there is a N_tr x d_i feature data
%test is a mx1 cell in each of which there is a N_te x d_i feature data
%trainlabel is a N_tr x 1 label vector for training set

if length(train)~=length(test)
    error('the number of types in training data is not same with test data');
end



for i = 1:length(train)
    if size(train{i},1)~=length(trainlabel)
        error('the number of points in %d -th data type is not same with the number of labels');
    end
    if size(train{i},2)~=size(test{i},2)
        error('The %d-th data has different dimensions in train and test', i);
    end
    Data = Standard_Normalization([train{i};test{i}]);
    Distance = dist2(Data,Data);
    Wall{i} = affinityMatrix(Distance, K, alpha);
end

%next, we fuse all the graphs
% then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(Wall, K, T);

N = size(W,1);Nl = size(train{1},1);
U = unique(trainlabel);
C = length(U);%%number of classes
Alllabel = zeros(N,C);%% 1-of-C labels
for i = 1: C
    Alllabel(trainlabel==U(i),i) = 1;
end


Alllabel = simpleLP(W,Alllabel, 1000);

[a,b] = max(Alllabel(1+Nl:end,:), [], 2);

Predicted_labels = U(b);




end

