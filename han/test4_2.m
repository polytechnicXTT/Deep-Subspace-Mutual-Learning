% 使用神经网络进行数据融合
clc
clear
close all

model = "model7";
data_view = "1";
view = "view_2";

% model_dir = 'D:\train\xlsx\OVARY\model10\';
model_dir = 'D:\PyCharm\PyCharm_Projects\MvSCN\xlsx\COAD\'+model+'\';
epoch = 20;
C = 3;
save_dir =  63;
K =            30;        %number of neighbors, usually (10~30)
alpha =     0.5;  %hyperparameter, usually (0.3~0.8)
s_c = [];
surv=importdata('D:\PyCharm\PyCharm_Projects\MvSCN\data\COAD\COLON_Survival.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;

    disp("当前运行的为：alpha:"+alpha+",K:"+K)
for i = 0:epoch-1
    Data = xlsread(model_dir + "W"+ i + "_" + data_view + ".xlsx", 'B2:JR278');
    Data = Standard_Normalization(Data);
    Dist = dist2(Data,Data);
    W = affinityMatrix(Dist, K, alpha);
    group = SpectralClustering(W,C);
    silhouette1=New_silhouette(W,group);
    s_value = "epoch "+ i + " : "+silhouette1;
    s_c(i+1) = silhouette1;
    %disp(s_value);
    outpart1="D:/PyCharm/PyCharm_Projects/MvSCN/han/test_results/COAD/single_view_res/"+view+"/"+model+"/"+model+"_param_comb/"+ save_dir +"/txt/colon_test"+i;
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
max_silhouette = max(s_c);
% disp(max_silhouette)
%disp("done!")

