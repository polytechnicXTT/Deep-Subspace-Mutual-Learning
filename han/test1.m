% 对模型生成的C1,C2,C3执行SNF融合
% 1. 与test2中的模型生成的C_fusion的聚类结果做对比实验
% 2. 与SNF方法本身的结果做对比实验

clc
clear
close all

% model_dir = 'D:\train\xlsx\OVARY\model10\';
model_dir = 'D:\PyCharm\Pycharm Projects\MvSCN\xlsx\KRCCC\model20\';
epoch = 20;
C = 3;
K =        14 ;%number of neighbors, usually (10~30)
alpha = 0.8 ;%hyperparameter, usually (0.3~0.8)
T = 20;%Number of Iterations, usually (10~20)

s_c = [];

surv=importdata('D:\PyCharm\Pycharm Projects\MvSCN\data\KRCCC\KIDNEY_Survival.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;

for i = 0:epoch-1
    Data1 = xlsread(model_dir + "W"+ i + "_0.xlsx",'B2:JR278');
    Data2 = xlsread(model_dir + "W"+ i + "_1.xlsx",'B2:JR278');
    Data3 = xlsread(model_dir + "W"+ i + "_2.xlsx",'B2:JR278');

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
    W = SNF({W1,W2,W3}, K, T);
    save("W.mat", "W")

    group = SpectralClustering(W,C);
    silhouette1=New_silhouette(W,group);
    s_value = "epoch "+ i + " : "+silhouette1;
    s_c(i+1) = silhouette1;
    disp(s_value);
    outpart1="D:/PyCharm/Pycharm Projects/MvSCN/han/test_results/KRCCC/test1/kidney_test"+i;
    outpart2=".txt";
    r_name=sprintf('%s%s',outpart1,outpart2);
    fid=fopen(r_name,'w');
    fprintf(fid,'patient_id     survival     death       Labels');
    fprintf(fid,'\n');
    for i=1:1:size(case_surv,1)
      fprintf(fid,'%s	%d	%d	%d',case_surv{i,1},data_surv(i,2),data_surv(i,1),group(i));
      fprintf(fid,'\n');
    end
    fclose(fid);
    disp(outpart1+outpart2+" saved")
end
max_silhouette = max(s_c)
% disp(max_silhouette)
disp("done!")

