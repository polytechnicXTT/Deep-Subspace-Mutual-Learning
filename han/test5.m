%¶ÔC1/C2/C3Ö´ÐÐÆ×¾ÛÀà
clc
clear
close all

% model_dir = 'D:\train\xlsx\OVARY\model_7\';
model_dir = 'D:\PyCharm\PyCharm Projects\MvSCN\xlsx\OVARY\Contrast\model1\';
epoch = 20;
C = 4;
K =          23 ;                  %number of neighbors, usually (10~30)
alpha =   0.5 ;             %hyperparameter, usually (0.3~0.8)

s_c = [];

surv=importdata('D:\PyCharm\PyCharm Projects\MvSCN\data\OVARY\OVARY_Survival.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;

for i = 0:epoch-1
    Data1 = xlsread(model_dir + "W"+ i + "_0.xlsx",'B2:JR278');
    Data1 = Standard_Normalization(Data1);
    %%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
    %if the data is discrete, we recommend the users to use chi-square distance
    Dist1 = dist2(Data1,Data1);
    %%%next, construct similarity graphs
    W1 = affinityMatrix(Dist1, K, alpha);

    group = SpectralClustering(W1,C);
    silhouette1=New_silhouette(W1,group);
    s_value = "epoch "+ i + " : "+silhouette1;
    s_c(i+1) = silhouette1;
    disp(s_value);
    outpart1="D:\PyCharm\PyCharm Projects\MvSCN\han\test_results\OVARY\test5\ovary_test"+i;
    outpart2=".txt";
    r_name=sprintf('%s%s',outpart1,outpart2);
    fid=fopen(r_name,'w');
    fprintf(fid,'patient_id     death     survival       Labels');
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

