% 使用神经网络进行数据融合
clc
clear
close all

cancer_type = "GBM";
save_file = "glio";
model = "model";
%data_view = "2";
data_view = "GLIO_Mirna_Expression";
view = "view_3";

% model_dir = 'D:\train\xlsx\OVARY\model10\';
% model_dir = 'D:\PyCharm Project\MvSCN\xlsx\COAD\'+model+'\';
data_dir = 'D:\PyCharm\Pycharm_Projects\DSML\data\'+cancer_type+'\';
epoch = 20;
C = 3;
K =        10;        %number of neighbors, usually (10~30)
alpha = 0.3;  %hyperparameter, usually (0.3~0.8)
%s_c = []
save_dir = 1;
surv=importdata('D:\PyCharm\PyCharm_Projects\DSML\data\GBM\GLIO_Survival.xlsx');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
%surv.data = surv.data(:,[2,1]);
surv.data = surv.data();
data_surv=surv.data;
for apha = 1:6
for K =10:1:30
    disp("当前运行的为：alpha:"+alpha+",K:"+K);

    % Data = xlsread(model_dir + "W"+ i + "_"+data_view+".xlsx",'B2:JR278');
    Data = xlsread(data_dir + data_view + ".xlsx", 'B2:XA625');
    Data = Standard_Normalization(Data);
    Dist = dist2(Data,Data);
    W = affinityMatrix(Dist, K, alpha);
    group = SpectralClustering(W,C);

    outpart1="D:/PyCharm/PyCharm_Projects/DSML/cluster_result/"+cancer_type+"/spectral_cluster_res/"+view+"/"+model+"/"+model+"_param_comb/"+ save_dir +"/txt/"+save_file+"_test.txt";
    r_name=sprintf('%s',outpart1);
    fid=fopen(r_name,'w');
    fprintf(fid,'patient_id     survival     death       Labels');
    fprintf(fid,'\n');
    for i=1:1:size(case_surv,1)
      fprintf(fid,'%s	%d	%d	%d',case_surv{i,1},data_surv(i,1),data_surv(i,2), group(i));
      fprintf(fid,'\n');
    end
    fclose(fid);
    disp(outpart1+" saved");

% max_silhouette = max(s_c);
% disp(max_silhouette)
% disp("done!")
save_dir = save_dir + 1;
end
alpha = alpha + 0.1;
end
