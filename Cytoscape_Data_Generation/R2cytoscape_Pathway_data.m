%% create cytoscape compatible edge and node attribute tables
clc
clear

D = dir('AROM_Paper_Codes\Cytoscape_Data_Generation\Sample_Data_for_pathway_analysis\*.xlsx');
if isempty(D)
    error('No xlsx files to process');
end

for k=1:length(D)
    %% Handle temporary file;just skip them
    if contains(D(k).name, '~')
        continue
    end
    disp(strcat('processing for file = ',D(k).name));
    T = readtable(strcat('AROM_Paper_Codes\Cytoscape_Data_Generation\Sample_Data_for_pathway_analysis\', D(k).name));
    d = T(:,[2, 6, 8, 9]);
    cyto = [];
    %% create edge attribute table
    for i=1:height(d)-1
        temp = strsplit(d.geneID{i}, '/');
        for j=i+1:height(d)
            temp2 = strsplit(d.geneID{j}, '/');
            %         fprintf('Common genes between %s and %s\n', d.Description{i}, d.Description{j});
            common_genes = intersect(temp, temp2);
            if ~isempty(common_genes)
                cyto = [cyto; [d.Description(i)  d.Description(j) num2cell(length(common_genes))]];
            end
        end
    end    
    cyto = cell2table(cyto, 'VariableNames', {'Source', 'Target', 'Edge_Attr'});
    
    %% Map edge attributes to the closed interval [desired_lower_bound, desired_upper_bound]
    cyto_node_attr = [];
    current_lower_bound = 1;
    current_upper_bound = max(cyto.Edge_Attr);
    desired_lower_bound = current_lower_bound;
    desired_upper_bound = 4;
    Mapped_Edge_Attr = map_numbers_to_an_interval(cyto.Edge_Attr, current_lower_bound, current_upper_bound, desired_lower_bound, desired_upper_bound);
    cyto = addvars(cyto, Mapped_Edge_Attr, 'After', 'Edge_Attr')
    
    %% create node attribute table
    cyto_node_attr = d(:,[1,2,end]);
    current_lower_bound = min(cyto_node_attr.Count);
    current_upper_bound = max(cyto_node_attr.Count);
    desired_lower_bound = 15;
    desired_upper_bound = 25;
    Mapped_Node_Attr = map_numbers_to_an_interval(cyto_node_attr.Count, current_lower_bound, current_upper_bound, desired_lower_bound, desired_upper_bound);
    cyto_node_attr = addvars(cyto_node_attr, Mapped_Node_Attr, 'After', 'Count');
    cyto_node_attr
    [~,ia,~] = unique(cyto_node_attr.Count);
    cyto_node_attr(ia,:)
    writetable(cyto, strcat('AROM_Paper_Codes\Cytoscape_Data_Generation\Cytoscape_NW_Tab_Pathway\Cyto_Edge_', regexprep(D(k).name, 'xlsx', 'txt')), 'Delimiter', '\t', 'WriteVariableNames', true);
    writetable(cyto_node_attr, strcat('AROM_Paper_Codes\Cytoscape_Data_Generation\Cytoscape_NW_Tab_Pathway\Cyto_Node_', regexprep(D(k).name, 'xlsx', 'txt')), 'Delimiter', '\t', 'WriteVariableNames', true);
end
