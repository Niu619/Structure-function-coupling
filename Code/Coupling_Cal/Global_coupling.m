%% Calcualte whole-brain coupling
sub = dir('F:\Cam_CAN\SC-FC\structure\GM_vol');
sub(1:2) = [];    
gray_mask = spm_vol('D:\Cam_CAN_fMRI\mask\GreyMask_vsize_1.5.nii');      % 给每个人灰质体积图谱卡上灰质mask
[X,~] = spm_read_vols(gray_mask);
for i = 1:length(sub)
    struc = spm_vol(strcat('F:\Cam_CAN\SC-FC\structure\GM_vol\',sub(i).name,'\','smwp1',sub(i).name,'_T1w.nii'));
    [X_struc,~] = spm_read_vols(struc);
    X_struc = X.*X_struc;
    spm_write_vol(struc,X_struc);
    clear struc X_struc
end

clear 
sub = dir('F:\Cam_CAN\SC-FC\structure\GM_vol');
sub(1:2) = [];    
total_coup = zeros(length(sub),1);
for i = 1:length(sub)
    struc = spm_vol(strcat('F:\Cam_CAN\SC-FC\structure\GM_vol\',sub(i).name,'\','smwp1',sub(i).name,'_T1w.nii'));
    [X_struc,~] = spm_read_vols(struc);
    func = spm_vol(strcat('F:\Cam_CAN\SC-FC\function\ALFF\',sub(i).name,'\','ALFFMap_',sub(i).name,'.nii'));
    [X_func,~] = spm_read_vols(func);
    X_struc = reshape(X_struc,121*145*121,1);
    X_func = reshape(X_func,121*145*121,1);
    [PDF_struc,~] = gretna_PDF(X_struc, 2^7, 'ksdensity');
    [PDF_func,~] = gretna_PDF(X_func, 2^7, 'ksdensity');
    total_coup(i,1) = gretna_KLDs(PDF_struc,PDF_func);
    clear struc X_struc func X_func PDF_struc PDF_func
end