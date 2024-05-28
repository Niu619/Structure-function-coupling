%% Split-half validation analysis-----------------------------------------------------------
% age
coupling_data = xlsread('F:\Cam_CAN\SC-FC\coupling\coupling246.xlsx','Sheet1','c2:in423'); 
x2 = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AN2:AN212');
coupling = coupling_data(x2(:),:);    
covariate(:,1) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AD2:AD212');
covariate(:,2) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AA2:AA212');
covariate(:,3) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AB2:AB212');
covariate(:,4:8) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AE2:AI212');
[beta,ResMS,T,P] = dong_multi_regress(coupling,covariate);
beta_age = beta(2,:)';
beta_ageage = beta(3,:)';
p_age = P(2,:)';
p_ageage = P(3,:)';
t_age = T(2,:)';
t_ageage = T(3,:)';
clear beta P T ResMS
[pID_age,~] = FDR(p_age,0.05);
a_age = find(p_age < pID_age);
b_age = find(p_age >= pID_age);
t_age(b_age(:)) = 0;
[pID_ageage,~] = FDR(p_ageage,0.05);
a_ageage = find(p_ageage < pID_ageage);
b_ageage = find(p_ageage >= pID_ageage);
t_ageage(b_ageage(:)) = 0;

% partial correlation fluid intelligence
coupling_age = coupling(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AE2:AI212');
fluid = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','X2:X212');
[r,p] = partialcorr(coupling_age,fluid,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(1,a())';
r(a(),1);

% mediation analysis
age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','Y2:Y212');
fluid = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','W2:W212');
[~,res_age,~,~] = dong_regress(age,covariate_age);
[~,res_fluid,~,~] = dong_regress(fluid,covariate_age);
norm_age = zscore(res_age);
norm_fluid = zscore(res_fluid);
result_beta = zeros(16,5);
result_p = zeros(16,5);
for i = 1:16
    data = coupling(:,i);
    [~,res_data,~,~] = dong_regress(data,covariate_age);
    norm_data = zscore(res_data);
    [~, stats] = mediation(norm_age, norm_fluid, norm_data,'boottop','plots','doCIs', 'bootsamples', 10000);
    result_beta(i,:) = stats.mean;
    result_p(i,:) = stats.p;
    clear data res_data norm_data stats
end

% partial correlation reaction time
coupling_data = xlsread('F:\Cam_CAN\SC-FC\coupling\coupling246.xlsx','Sheet1','c2:in423'); 
x2 = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AN2:AN212');
coupling = coupling_data(x2(:),:);    
mRT = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AM2:AM212');
site_RT = find(mRT==0);
mRT(site_RT(),:) = [];
mRT = zscore(mRT);
coupling(site_RT(),:) = [];
a_age = xlsread('C:\Users\dell\Desktop\subregion_network_Yeo.xlsx','Sheet1','G90:G105'); 
coupling_age = coupling(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','S-Fvalidation','AE2:AI212');
covariate_age(site_RT(),:) = [];
[r,p] = partialcorr(coupling_age,mRT,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(a(),1);
r(a(),1);

% mediation analysis reaction time
age = xlsread('C:\Users\dell\Desktop\行为数据.xlsx','S-Fvalidation','Y2:Y212');
age(site_RT()) = [];
mRT = xlsread('C:\Users\dell\Desktop\行为数据.xlsx','S-Fvalidation','AM2:AM212');
mRT(site_RT()) = [];
[~,res_age,~,~] = dong_regress(age,covariate_age);
[~,res_mRT,~,~] = dong_regress(mRT,covariate_age);
norm_age = zscore(res_age);
norm_mRT = zscore(res_mRT);
result_beta = zeros(16,5);
result_p = zeros(16,5);
for i = 1:16
    data = coupling(:,i);
    [~,res_data,~,~] = dong_regress(data,covariate_age);
    norm_data = zscore(res_data);
    [~, stats] = mediation(norm_age, norm_mRT, norm_data,'boottop','plots','doCIs', 'bootsamples', 10000);
    result_beta(i,:) = stats.mean;
    result_p(i,:) = stats.p;
    clear data res_data norm_data stats
end
