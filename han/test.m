clc
clear
close all

K = 20;%number of neighbors, usually (10~30)
alpha = 0.4; %hyperparameter, usually (0.3~0.8)
T = 20; %Number of Iterations, usually (10~20)

model_dir = 'E:\PyCharm Project\train\xlsx\OVARY\model10\';
% model_dir = 'D:\test\xlsx\model10\';
epoch = 4;
Data1 = xlsread(model_dir + "W"+ epoch + "_0.xlsx",'B2:JR278');
Data2 = xlsread(model_dir + "W"+ epoch + "_1.xlsx",'B2:JR278');
Data3 = xlsread(model_dir + "W"+ epoch + "_2.xlsx",'B2:JR278');

% Data1 = xlsread("D:\train\best\ovgay\ovgay_ge.xlsx",'B2:JR7499')';
% Data2 = xlsread("D:\train\best\ovgay\ovgay_me.xlsx",'B2:JR22287')';
% Data3 = xlsread("D:\train\best\ovgay\ovgay_mi.xlsx",'B2:JR340')';

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

% clear all;
% clc;

surv=importdata('D:/data/ovarian_cancer_data/survival_data.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;
% load('W');
C = 4;
group = SpectralClustering(W,C);
silhouette1=New_silhouette(W,group);
displayClusters(W,group);
%title('New method');
outpart1='./test_results/ovary_test';
outpart2='.txt';
r_name=sprintf('%s%s',outpart1,outpart2);
fid=fopen(r_name,'w');
fprintf(fid,'patient_id	death	suvival	Labels');
fprintf(fid,'\n');
for i=1:1:size(case_surv,1)
  fprintf(fid,'%s	%d	%d	%d',case_surv{i,1},data_surv(i,2),data_surv(i,1),group(i));
  fprintf(fid,'\n');
end
fclose(fid);
