% 对原始数据单独执行谱聚类
% 与test4中对应的C1/C2/C3单独进行SNF再进行聚类结果做对比实验
clc
clear
close all

% model_dir = 'D:\train\xlsx\OVARY\model10\';
%model_dir = 'D:\Pycharm Project\MvSCN\xlsx\BREAST\model\';
data_dir = 'D:\PyCharm\PyCharm Projects\MvSCN\data\KRCCC\';
epoch = 20;
C = 3;
K = 11;%number of neighbors, usually (10~30)
alpha = 0.3;%hyperparameter, usually (0.3~0.8)

surv=importdata('D:\PyCharm\PyCharm Projects\MvSCN\data\KRCCC\KIDNEY_Survival.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;

Data1 = xlsread(data_dir + "Kidney_methy.xlsx", 'B2:JR278');
Data1 = Standard_Normalization(Data1);
    
%%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
%if the data is discrete, we recommend the users to use chi-square distance
Dist1 = dist2(Data1,Data1);
%%%next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha);

group = SpectralClustering(W1,C);
silhouette1=New_silhouette(W1,group);
s_value = "S_value : "+silhouette1;
disp(s_value);
outpart1="D:/PyCharm/PyCharm Projects/MvSCN/han/test_results/KRCCC/test3/kidney_test.txt";
r_name=sprintf('%s',outpart1);
fid=fopen(r_name,'w');
fprintf(fid,'patient_id     survival     death       Labels');
fprintf(fid,'\n');
for i=1:1:size(case_surv,1)
  fprintf(fid,'%s	%d	%d	%d',case_surv{i,1},data_surv(i,2),data_surv(i,1),group(i));
  fprintf(fid,'\n');
end
fclose(fid);
disp(outpart1+" saved")

%max_silhouette = max(s_c)
% disp(max_silhouette)
disp("done!")

