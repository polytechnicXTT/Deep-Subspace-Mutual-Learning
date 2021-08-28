 %¶ÔC1/C2/C3Ö´ÐÐÆ×¾ÛÀà
clc
clear
close all

cancer_type = "COAD";
C = 3;
model = "model12";
data_view = "1";
view = "view_2";
save_file_name = "colon";


% model_dir = 'D:\train\xlsx\OVARY\model_7\';
model_dir = 'D:\PyCharm\Pycharm_Projects\DSML\model\'+cancer_type+'\'+model+'\';
epoch = 20;
save_dir = 63;
K =           30 ;%number of neighbors, usually (10~30)
alpha = 0.5 ;%hyperparameter, usually (0.3~0.8)

% s_c = [];

surv=importdata('D:\PyCharm\Pycharm_Projects\DSML\data\COAD\COLON_Survival.xlsx');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;

for i = 0:epoch-1
    Data1 = xlsread(model_dir + "W"+ i + "_"+data_view+".xlsx",'B2:XA625');
    Data1 = Standard_Normalization(Data1);
    %%%Calculate the pair-wise distance; If the data is continuous, we recommend to use the function "dist2" as follows;
    %if the data is discrete, we recommend the users to use chi-square distance
    Dist1 = dist2(Data1,Data1);
    %%%next, construct similarity graphs
    W1 = affinityMatrix(Dist1, K, alpha);

    group = SpectralClustering(W1,C);
%     silhouette1=New_silhouette(W1,group);
%     s_value = "epoch "+ i + " : "+silhouette1;
%     s_c(i+1) = silhouette1;
%     disp(s_value);
    outpart1="D:/PyCharm/Pycharm_Projects/DSML/cluster_result/"+cancer_type+"/single_view_res/"+view+"/"+model+"/"+model+"_param_comb/"+ save_dir +"/txt/"+save_file_name+"_test"+i;
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
% max_silhouette = max(s_c)
% disp(max_silhouette)
disp("done!")

