%% This code integrates MALDI imaging and LC-MS datasets to extract meaningful information
%% This code first fixes the Scills Lab masslists corresponding to peptide clusters using entire IMS spectrum
%% Then searches those peptides in peptide LIMMA(done on LC-MS data) file
%% After that it resolves the ambiguity: "one IMS peptide could match multple LC-MS peptides within the accuracy(± 0.1 Da) of IMS measurements" using MLP scoring


clc
clear
close all
format longg

%% Find how many excel files are in the specific folder
dirname = 'Sample_Data_For_MLP';
D = dir(strcat(dirname, '/*.xlsx'));
if isempty(D)
    error('No xlsx files to process');
end

%% Create folder for output
if ~exist(strcat(dirname,'/mass_lists_after_Deisotoping'), 'dir')
    mkdir (strcat(dirname,'/mass_lists_after_Deisotoping'));
end

if ~exist(strcat(dirname,'/Output'), 'dir')
    mkdir (strcat(dirname,'/Output'));
end

Log_file = strcat(dirname,'/Output/FAILURE_LOG_FILE.txt');
fileID = fopen(Log_file,'a');
tol = 0.1;
for i=1:length(D)
    
    %% Handle temporary file;just skip them
    if contains(D(i).name, '~')
        continue
    end
    filename_with_path = strcat(dirname, '/', D(i).name);
    [IMS_mass, ~, ~] = xlsread(filename_with_path);
    temp = IMS_mass(:,1);
    
    %% deisotoping of IMS_mass: for a mass 'm', remove all masses that fall in the interval, [(m+0.85),(m+1.15)]
    IMS_mass = MLP_Helper_Functions.deisotoping_IMS(temp);
    IMS_mass_bak = IMS_mass;
    
    %% Second Round of Deisotoping to fix Scills Lab's wrong Peak assignment
    Spectrum = table2cell(readtable('WT_R2_WMS_cleaned.csv', 'Delimiter', ','));
    Spectrum = cell2mat(Spectrum);
    IMS_round = round(IMS_mass);
    cm = MLP_Helper_Functions.find_consecutive_numbers(IMS_round);
    if ~isempty(cm)
        mass_corrected = MLP_Helper_Functions.replace_with_masses_correspoding_to_correct_peak(cm, Spectrum);
        [~, idx] = ismember(cm, IMS_round);
        IMS_mass(idx) = mass_corrected;
        IMS_mass = MLP_Helper_Functions.deisotoping_IMS(IMS_mass);
        file_name = strsplit(D(i).name, '.xlsx');
        file_name = strcat(dirname,'/mass_lists_after_Deisotoping/', file_name{1}, '.xls');
        xlswrite(file_name, IMS_mass);
    else
        file_name = strsplit(D(i).name, '.xlsx');
        file_name = strcat(dirname,'/mass_lists_after_Deisotoping/', file_name{1}, '.xls');
        xlswrite(file_name, IMS_mass);
    end
    
    %% Check if the IMS_mass to be searched in AROM or WT LIMMA data
    if contains(D(i).name, 'AROM')
        
        %% LogFC > 0 AND Pvalue DON'T CARE
        MLP_Helper_Functions.serach_in_LIMMA_results( strcat(dirname, '/Results_LogFC_1_point_5_Peptide_File/final_data_treatment_FC_grt_0_Pval_DC.txt'), ...
            '_LIMMA_FC_grt_0_Pvalue_DC.xls', ...
            'FC_grt_0_Pvalue_DC:', ...
            fileID, ...
            filename_with_path, ...
            dirname, ...
            tol, ...
            IMS_mass);
        
        %% LogFC > 1.5 AND Pvalue < 0.05
        MLP_Helper_Functions.serach_in_LIMMA_results(strcat(dirname, '/Results_LogFC_1_point_5_Peptide_File/final_data_treatment.txt'), ...
            '_LIMMA_FC_grt_1p5_Pvalue_less_0p05.xls', ...
            'FC_grt_1p5_Pvalue_less_0p05:', ...
            fileID, ...
            filename_with_path, ...
            dirname, ...
            tol, ...
            IMS_mass);
        
        
    else
        %% LogFC < 0 AND Pvalue DON'T CARE
        MLP_Helper_Functions.serach_in_LIMMA_results(strcat(dirname, '/Results_LogFC_1_point_5_Peptide_File/final_data_control_FC_less_0_Pval_DC.txt'), ...
            '_LIMMA_FC_less_0_Pvalue_DC.xls', ...
            'FC_less_0_Pvalue_DC:', ...
            fileID, ...
            filename_with_path, ...
            dirname, ...
            tol, ...
            IMS_mass);
        
        
        %% LogFC < 1.5 AND Pvalue < 0.05
        MLP_Helper_Functions.serach_in_LIMMA_results(strcat(dirname, '/Results_LogFC_1_point_5_Peptide_File/final_data_control.txt'), ...
            '_LIMMA_FC_less_1p5_Pvalue_less_0p05.xls', ...
            'FC_less_1p5_Pvalue_less_0p05:', ...
            fileID, ...
            filename_with_path, ...
            dirname, ...
            tol, ...
            IMS_mass);
        
        
    end
end



%% Resolves the ambiguity: "one IMS peptide could match multple LC-MS peptides within the accuracy(Â± 0.1 Da) of IMS measurements" using MLP scoring
D = dir(strcat(dirname, '/Output/*.xls'));

if isempty(D)
    error('No xls files to process');
end

for i=1:length(D)
    %% Handle temporary file;just skip them
    if contains(D(i).name, '~')
        continue
    end
    filename_with_path = strcat(dirname, '/Output/', D(i).name);
    [~,~,T_bak] = xlsread(filename_with_path);
    T = T_bak(2:end,:);
    tol = 0.2; %% increased tol to make sure boundary is achieved
    k = 1;
    Top_hits_with_uniq = [];
    second_hits_with_uniq = [];
    Mass = cell2mat(T(:,1));
    
    %% MLP scoring
    while k <= length(Mass)
        m = Mass(k);
        idx = find((Mass >= (m-tol)) & (Mass <= (m+tol)));
        grouped_data = T(idx,:);
        MLP = (str2double(T(idx, end-3)).*abs(str2double(T(idx, end-1))))./str2double(T(idx, end));
        max_MLP = max(MLP);
        MLP = MLP/max_MLP(1);
        k = k + length(idx);
        [MLP_sorted, I] = sort(MLP, 'descend');
        grouped_data = [grouped_data(I,:) num2cell(MLP_sorted)];
        if(size(grouped_data,1) == 1)
            Top_hits_with_uniq = [Top_hits_with_uniq; grouped_data];
        else
            Top_hits_with_uniq = [Top_hits_with_uniq; grouped_data(1,:)];
            second_hits_with_uniq = [second_hits_with_uniq; grouped_data(2,:)];
        end
    end
    
    %% Write Top Hits
    if ~exist(strcat(dirname, '/Output/MLP'), 'dir')
        mkdir (strcat(dirname, '/Output/MLP'));
    end
    [~,filename_without_ext,~] = fileparts(filename_with_path);
    outfileTop = strcat(dirname, '/Output/MLP/', 'Top_MLP_hits_', filename_without_ext, '.xls');
    xlswrite(outfileTop, [[T_bak(1,:) {'MLP Score'}]; Top_hits_with_uniq]);
    
    %% Write Second Hits
    if ~isempty(second_hits_with_uniq)
        outfileSecond = strcat(dirname, '/Output/MLP/', 'Second_MLP_hits_', filename_without_ext, '.xls');
        xlswrite(outfileSecond, [[T_bak(1,:) {'MLP Score'}]; second_hits_with_uniq]);
    end
end

fclose('all');