clear all;
clc;

%surv=importdata('E:\PyCharm Project\train\data\ovarian cancer data\survival data.txt');
surv=importdata('D:\PyCharm Projects\MvSCN\data\COAD\COLON_Survival.txt');
STR_surv = surv.textdata(2:size(surv.textdata,1),1);
case_surv=arrayfun(@(k) STR_surv{k}(1:12),1:length(STR_surv),'UniformOutput',0)';
surv.data = surv.data(:,[2,1]);
data_surv=surv.data;
load('W');
group = SpectralClustering(W,6);
silhouette1=New_silhouette(W,group);
displayClusters(W,group);
%title('New method');
outpart1='D:/PyCharm Projects/MvSCN/han/test_results/colon_test';
outpart2='.txt';
r_name=sprintf('%s%s',outpart1,outpart2);
fid=fopen(r_name,'w');
fprintf(fid,'patient_id     survival    death     Labels');
fprintf(fid,'\n');
for i=1:1:size(case_surv,1)
  fprintf(fid,'%s	%d	%d	%d',case_surv{i,1},data_surv(i,2),data_surv(i,1),group(i));
  fprintf(fid,'\n');
end
fclose(fid);
