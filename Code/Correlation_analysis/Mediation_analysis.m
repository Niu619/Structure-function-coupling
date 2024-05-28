%% Mediation analysis
%% BN ATLAS
% fluid intelligence
SFG_R_7_3 = coupling_index(:,6);
STG_L_6_3 = coupling_index(:,73);
STG_R_6_3 = coupling_index(:,74);
STG_R_6_6 = coupling_index(:,80);
ITG_R_7_2 = coupling_index(:,92);
PhG_L_6_5 = coupling_index(:,117);
SPL_L_5_4 = coupling_index(:,131);
IPL_L_6_2 = coupling_index(:,137);
INS_R_6_1 = coupling_index(:,164);
INS_L_6_3 = coupling_index(:,167);
LOcC_L_4_2 = coupling_index(:,201);
BG_L_6_1 = coupling_index(:,219);
BG_L_6_3 = coupling_index(:,223);
Tha_R_8_1 = coupling_index(:,232);
Tha_L_8_6 = coupling_index(:,241);
Tha_R_8_6 = coupling_index(:,242);

age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','D2:D423');
FI = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','B2:B423');
covariate = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
[~,res_age,~,~] = dong_regress(age,covariate);
[~,res_FI,~,~] = dong_regress(FI,covariate);
[~,res_SFG_R_7_3,~,~] = dong_regress(SFG_R_7_3,covariate);
[~,res_STG_L_6_3,~,~] = dong_regress(STG_L_6_3,covariate);
[~,res_STG_R_6_3,~,~] = dong_regress(STG_R_6_3,covariate);
[~,res_STG_R_6_6,~,~] = dong_regress(STG_R_6_6,covariate);
[~,res_ITG_R_7_2,~,~] = dong_regress(ITG_R_7_2,covariate);
[~,res_PhG_L_6_5,~,~] = dong_regress(PhG_L_6_5,covariate);
[~,res_SPL_L_5_4,~,~] = dong_regress(SPL_L_5_4,covariate);
[~,res_IPL_L_6_2,~,~] = dong_regress(IPL_L_6_2,covariate);
[~,res_INS_R_6_1,~,~] = dong_regress(INS_R_6_1,covariate);
[~,res_INS_L_6_3,~,~] = dong_regress(INS_L_6_3,covariate);
[~,res_LOcC_L_4_2,~,~] = dong_regress(LOcC_L_4_2,covariate);
[~,res_BG_L_6_1,~,~] = dong_regress(BG_L_6_1,covariate);
[~,res_BG_L_6_3,~,~] = dong_regress(BG_L_6_3,covariate);
[~,res_Tha_R_8_1,~,~] = dong_regress(Tha_R_8_1,covariate);
[~,res_Tha_L_8_6,~,~] = dong_regress(Tha_L_8_6,covariate);
[~,res_Tha_R_8_6,~,~] = dong_regress(Tha_R_8_6,covariate);
norm_age = zscore(res_age);
norm_FI = zscore(res_FI);
norm_SFG_R_7_3 = zscore(res_SFG_R_7_3);
norm_STG_L_6_3 = zscore(res_STG_L_6_3);
norm_STG_R_6_3 = zscore(res_STG_R_6_3);
norm_STG_R_6_6 = zscore(res_STG_R_6_6);
norm_ITG_R_7_2 = zscore(res_ITG_R_7_2);
norm_PhG_L_6_5 = zscore(res_PhG_L_6_5);
norm_SPL_L_5_4 = zscore(res_SPL_L_5_4);
norm_IPL_L_6_2 = zscore(res_IPL_L_6_2);
norm_INS_R_6_1 = zscore(res_INS_R_6_1);
norm_INS_L_6_3 = zscore(res_INS_L_6_3);
norm_LOcC_L_4_2 = zscore(res_LOcC_L_4_2);
norm_BG_L_6_1 = zscore(res_BG_L_6_1);
norm_BG_L_6_3 = zscore(res_BG_L_6_3);
norm_Tha_R_8_1 = zscore(res_Tha_R_8_1);
norm_Tha_L_8_6 = zscore(res_Tha_L_8_6);
norm_Tha_R_8_6 = zscore(res_Tha_R_8_6);

[paths, stats] = mediation(norm_age, norm_FI, norm_Tha_R_8_6,'boottop','plots','doCIs', 'bootsamples', 10000);

% Reaction time
age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','D2:D423');
age(site_RT()) = [];
mRT = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','T2:T423');
mRT(site_RT()) = [];
[~,res_age,~,~] = dong_regress(age,covariate_age);
[~,res_mRT,~,~] = dong_regress(mRT,covariate_age);
norm_age = zscore(res_age);
norm_mRT = zscore(res_mRT);
result_beta = zeros(49,5);
result_p = zeros(49,5);
for i = 1:49
    data = coupling_index(:,i);
    [~,res_data,~,~] = dong_regress(data,covariate_age);
    norm_data = zscore(res_data);
    [~, stats] = mediation(norm_age, norm_mRT, norm_data,'boottop','plots','doCIs', 'bootsamples', 10000);
    result_beta(i,:) = stats.mean;
    result_p(i,:) = stats.p;
    clear data res_data norm_data stats
end

%% HOA ATLAS
% fluid intelligence
age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','D2:D423');
fluid = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','B2:B423');
[~,res_age,~,~] = dong_regress(age,covariate_age);
[~,res_fluid,~,~] = dong_regress(fluid,covariate_age);
norm_age = zscore(res_age);
norm_fluid = zscore(res_fluid);
result_beta = zeros(19,5);
result_p = zeros(19,5);
for i = 1:19
    data = coupling_index(:,i);
    [~,res_data,~,~] = dong_regress(data,covariate_age);
    norm_data = zscore(res_data);
    [~, stats] = mediation(norm_age, norm_fluid, norm_data,'boottop','plots','doCIs', 'bootsamples', 10000);
    result_beta(i,:) = stats.mean;
    result_p(i,:) = stats.p;
    clear data res_data norm_data stats
end

% reaction time
age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','D2:D423');
age(site_RT()) = [];
mRT = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','T2:T423');
mRT(site_RT()) = [];
[~,res_age,~,~] = dong_regress(age,covariate_age);
[~,res_mRT,~,~] = dong_regress(mRT,covariate_age);
norm_age = zscore(res_age);
norm_mRT = zscore(res_mRT);
result_beta = zeros(19,5);
result_p = zeros(19,5);
for i = 1:19
    data = coupling_index(:,i);
    [~,res_data,~,~] = dong_regress(data,covariate_age);
    norm_data = zscore(res_data);
    [~, stats] = mediation(norm_age, norm_mRT, norm_data,'boottop','plots','doCIs', 'bootsamples', 10000);
    result_beta(i,:) = stats.mean;
    result_p(i,:) = stats.p;
    clear data res_data norm_data stats
end
