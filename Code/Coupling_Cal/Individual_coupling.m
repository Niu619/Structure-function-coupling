%% Calculate individual coupling    
Data_path = 'F:\Cam_CAN\SC-FC\structure\data_path.txt';
File_filter = 'smwp1';
Template_path = 'F:\Cam_CAN\SC-FC\structure\Template_path.txt';
Template_filter = 'BN_Atlas';
RoiIndex = [1:246];
Output_path = 'F:\Cam_CAN\SC-FC\structure\vol_signal';
Level = 'all';
gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level);
clear
% 
gretna_gen_data_path({'F:\Cam_CAN\SC-FC\function\ALFF'},'F:\Cam_CAN\SC-FC\function','data_path.txt');
Data_path = 'F:\Cam_CAN\SC-FC\function\data_path.txt';
File_filter = 'ALFFMap';
Template_path = 'F:\Cam_CAN\SC-FC\function\Template_path.txt';
Template_filter = 'BN_Atlas';
RoiIndex = [1:246];
Output_path = 'F:\Cam_CAN\SC-FC\function\ALFF_signal';
Level = 'all';
gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level);

% 
struc = dir('F:\Cam_CAN\SC-FC\structure\vol_signal');
func = dir('F:\Cam_CAN\SC-FC\function\ALFF_signal');
struc(1:2) = []; func(1:2) = [];
coupling = cell(422,1);
for i = 1:422
    x = load(strcat('F:\Cam_CAN\SC-FC\structure\vol_signal\',struc(i).name));
    y = load(strcat('F:\Cam_CAN\SC-FC\function\ALFF_signal\',func(i).name));
    coupl = zeros(246,1);
    for j = 1:246
        [PDFx,~] = gretna_PDF(x.Sig_roi{1,j}, 2^7, 'ksdensity');
        [PDFy,~] = gretna_PDF(y.Sig_roi{1,j}, 2^7, 'ksdensity');
        coupl(j,1) = gretna_KLDs(PDFx,PDFy);
        clear PDFx PDFy
    end
    coupling{i} = coupl;        
    clear coupl x y 
end
total_coup = zeros(422,246);
for m = 1:422 
    total_coup(m,:) = coupling{m,1};
end   

%% HOA ATLAS

Data_path = 'F:\Cam_CAN\SC-FC\structure\data_path.txt';
File_filter = 'smwp1';
Template_path = 'F:\Cam_CAN\SC-FC\structure\Template_path.txt';
Template_filter = 'HOA';
RoiIndex = [1:110];
Output_path = 'F:\Cam_CAN\SC-FC\HOA_result\vol_signal';
Level = 'all';
gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level);
clear
% 
gretna_gen_data_path({'F:\Cam_CAN\SC-FC\HOA_result\ALFF'},'F:\Cam_CAN\SC-FC\HOA_result','data_path.txt');
Data_path = 'F:\Cam_CAN\SC-FC\HOA_result\data_path.txt';
File_filter = 'ALFFMap';
Template_path = 'F:\Cam_CAN\SC-FC\Template_path.txt';
Template_filter = 'HOA';
RoiIndex = [1:110];
Output_path = 'F:\Cam_CAN\SC-FC\HOA_result\ALFF_signal';
Level = 'all';
gretna_Vol_ROI_signal(Data_path, File_filter, Template_path, Template_filter, RoiIndex, Output_path, Level);
% 
struc = dir('F:\Cam_CAN\SC-FC\HOA_result\vol_signal');
func = dir('F:\Cam_CAN\SC-FC\HOA_result\ALFF_signal');
struc(1:2) = []; func(1:2) = [];

coupling = cell(422,1);
for i = 1:422
    x = load(strcat('F:\Cam_CAN\SC-FC\HOA_result\vol_signal\',struc(i).name));
    y = load(strcat('F:\Cam_CAN\SC-FC\HOA_result\ALFF_signal\',func(i).name));
    coupl = zeros(110,1);
    for j = 1:110
        [PDFx,~] = gretna_PDF(x.Sig_roi{1,j}, 2^7, 'ksdensity');
        [PDFy,~] = gretna_PDF(y.Sig_roi{1,j}, 2^7, 'ksdensity');
        coupl(j,1) = gretna_KLDs(PDFx,PDFy);
        clear PDFx PDFy
    end
    coupling{i} = coupl;
    clear coupl x y 
end
save('F:\Cam_CAN\SC-FC\HOA_result\coupling.mat','coupling');
