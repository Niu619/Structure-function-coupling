%% Relation to age and cognition
%% BN ATLAS
% regression analysis with age
coupling_index = xlsread('F:\Cam_CAN\SC-FC\coupling\coupling246.xlsx','Sheet1','c2:in423'); 
covariate(:,1) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','I2:I423');
covariate(:,2) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','F2:F423');
covariate(:,3) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','G2:G423');
covariate(:,4:8) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
[beta,ResMS,T,P] = dong_multi_regress(coupling_index,covariate);
beta_age = beta(2,:);
beta_ageage = beta(3,:);
p_age = P(2,:);
p_ageage = P(3,:);
t_age = T(2,:);
t_ageage = T(3,:);
clear beta P T ResMS
[pID_age,~] = FDR(p_age,0.05);
a_age = find(p_age < pID_age);
b_age = find(p_age >= pID_age);
t_age(b_age(:)) = 0;
[pID_ageage,~] = FDR(p_ageage,0.05);
a_ageage = find(p_ageage < pID_ageage);
b_ageage = find(p_ageage >= pID_ageage);
t_ageage(b_ageage(:)) = 0;

% partial correlation analysis with cognition
% fluid intelligence
coupling_age = coupling_index(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
fluid = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','C2:C423');
[r,p] = partialcorr(coupling_age,fluid,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(1,a())';
r(a(),1);

% reaction time
coupling_index = xlsread('F:\Cam_CAN\SC-FC\coupling\coupling246.xlsx','Sheet1','c2:in423');
mRT = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','T2:T423');
site_RT = find(mRT==0);
mRT(site_RT(),:) = [];
mRT = zscore(mRT);
coupling_index(site_RT(),:) = [];
a_age = xlsread('C:\Users\dell\Desktop\subregion_network_Yeo.xlsx','Sheet1','G2:G50'); 
coupling_age = coupling_index(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
covariate_age(site_RT(),:) = [];
[r,p] = partialcorr(coupling_age,mRT,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(a(),1);
r(a(),1);

%% HOA ATLAS
% age-------------------
coupling_index = xlsread('F:\Cam_CAN\SC-FC\HOA_result\coupling110.xlsx','Sheet1','c2:dh423');
covariate(:,1) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','I2:I423');
covariate(:,2) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','F2:F423');
covariate(:,3) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','G2:G423');
covariate(:,4:8) = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
[beta,ResMS,T,P] = dong_multi_regress(coupling_index,covariate);
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
coupling_age = coupling_index(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
fluid = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','C2:C423');
[r,p] = partialcorr(coupling_age,fluid,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(a(),1);
r(a(),1);

% partial correlation reaction time
coupling_index = xlsread('F:\Cam_CAN\SC-FC\HOA_result\coupling110.xlsx','Sheet1','c2:dh423');
mRT = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','T2:T423');
site_RT = find(mRT==0);
mRT(site_RT(),:) = [];
mRT = zscore(mRT);
coupling_index(site_RT(),:) = [];
a_age = xlsread('C:\Users\dell\Desktop\subregion_network_Yeo.xlsx','HOA','I8:I26'); 
coupling_age = coupling_index(:,a_age());
covariate_age = xlsread('C:\Users\dell\Desktop\behavior.xlsx','SC-FC','J2:N423');
covariate_age(site_RT(),:) = [];
[r,p] = partialcorr(coupling_age,mRT,covariate_age);
[pID,pN] = FDR(p,0.05);
a = find(p < pID);
b = find(p >= pID);
site = a_age(a(),1);
r(a(),1);


