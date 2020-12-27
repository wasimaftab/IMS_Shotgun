%% create cytoscape compatible edge table

clc
clear

D = dir('AROM_Paper_Codes\Cytoscape_Data_Generation\Sample_Data_for_GO_analysis\*.xlsx');
if isempty(D)
    error('No xlsx files to process');
end

for i=1:length(D)
    %% Handle temporary file;just skip them
    if contains(D(i).name, '~')
        continue
    end
    T = readtable(strcat('AROM_Paper_Codes\Cytoscape_Data_Generation\Sample_Data_for_GO_analysis\', D(i).name));
    
    d = T(:,[1, 2, 6, 8, 9]);
    cyto = [];
    base_size = 9;
    for k=1:height(d)
        temp = strsplit(d.geneID{k}, '/');
        for j = 1:length(temp)
            cyto = [cyto; [temp(j) d.Description(k) num2cell(d.Count(k))]];
        end
    end
    cyto = cell2table([cyto(:,1) cyto(:,2) cyto(:,end) ], 'VariableNames', {'Source' 'Target' 'Count'});
    writetable(cyto, strcat('AROM_Paper_Codes\Cytoscape_Data_Generation\Cytoscape_NW_Tab_GO\Cyto_', regexprep(D(i).name, 'xlsx', 'txt')), 'Delimiter', '\t', 'WriteVariableNames', true);    
end

% d = T(:,[1, 2, 6, 8, 9]);
% cyto = [];
% base_size = 9;
% for i=1:height(d)
%     temp = strsplit(d.geneID{i}, '/');
%     %     length(temp)
%     for j = 1:length(temp)
%         %         cyto = [cyto; [d.Description(i) temp(j) num2cell(length(temp))]];
%         cyto = [cyto; [temp(j) d.Description(i) num2cell(d.Count(i))]];
%     end
% %     keyboard
% end
%
% % cyto = cell2table([cyto(:,2) cyto(:,1) cyto(:,end) ], 'VariableNames', {'Source' 'Target' 'Count'});
% cyto = cell2table([cyto(:,1) cyto(:,2) cyto(:,end) ], 'VariableNames', {'Source' 'Target' 'Count'});
% % [~,I,~] = unique(cyto.Target, 'sorted');
% % temp = cyto(I,:)
% % [temp.Count map_numbers_to_an_interval(temp.Count, min(temp.Count), max(temp.Count), 15, 25)]
%
% writetable(cyto, 'AROM_Paper_Codes\Cytoscape_Data_Generation\Cytoscape_NW_Tab_GO\NW_Tab_Fig5B_GO.txt', 'Delimiter', '\t', 'WriteVariableNames', true);