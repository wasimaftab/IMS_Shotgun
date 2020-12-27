classdef MLP_Helper_Functions
    methods(Static)
        
        %% deisotoping of IMS_mass: for a mass 'm', remove all masses that fall in the interval, [(m+0.85),(m+1.15)]
        function IMS_mass = deisotoping_IMS(IMS_mass)
            idx_remove_IMS_mass = [];
            for i=1:length(IMS_mass)
                m = IMS_mass(i);
                lb = 0.85;
                ub = 1.15;
                idx_remove_IMS_mass = [idx_remove_IMS_mass; find((IMS_mass >= (m+lb)) & (IMS_mass <= (m+ub)))];
            end
            IMS_mass(idx_remove_IMS_mass) = [];
        end
        
        %% find consecutive numbers
        function d = find_consecutive_numbers(v)
            d = [];
            for i=1:length(v)-1
                if abs(v(i)-v(i+1)) == 1
                    d = [d; v(i:i+1)];
                end
            end
            
            d = unique(d);
        end
        
        %% This function returns the mass associated with the actual peak
        function  mass_corrected = replace_with_masses_correspoding_to_correct_peak(cm, Spectrum)
            mass_corrected = cm;
            for i = 1:length(cm)
                idx = find((Spectrum(:,1) >= cm(i)) & (Spectrum(:,1) <= cm(i)+1));
                idx2 = find(Spectrum(idx,2) == max(Spectrum(idx,2)));
                mass_corrected(i) = Spectrum(idx(idx2),1);
            end
            
        end
        
        
        %% Search mass with tolerance in Limma output obtained on Peptide file
        function Filtered_DS = search_mass_with_tol_in_Limma_Peptide_Data(tol, IMS_mass, Enriched_df_withhead)
            Enriched_df = Enriched_df_withhead(2:end,:);
            Pep_Mass = str2double(Enriched_df(:,1))+1.008;
            Enriched_df(:,1) = num2cell(Pep_Mass);
            % tol = 0.1;
            Filtered_DS = [];
            col_idx_mean_Inten = find(~cellfun(@isempty, regexpi(Enriched_df_withhead(1,:), 'mean')));
            for i=1:length(IMS_mass)
                idx_match = find((Pep_Mass>=IMS_mass(i)-tol) & (Pep_Mass<=IMS_mass(i)+tol));
                if ~isempty(idx_match)
                    temp = Enriched_df(idx_match,:);
                    [~,I] = sort(str2double(temp(:,col_idx_mean_Inten)), 'descend');
                    Filtered_DS = [Filtered_DS; temp(I,:)];
                end
            end
        end
        
        
        %% Search in LIMMA results
        function serach_in_LIMMA_results(LIMMA_file, out_file_end_part, Log_File_Tag, Log_File_Id, filename_with_path, dirname, tol, IMS_mass)
            LIMMA_pep = table2cell(readtable(LIMMA_file, 'Delimiter', '\t', 'ReadVariableNames', false));
            %% Search with tolerance
            Filtered_DS = MLP_Helper_Functions.search_mass_with_tol_in_Limma_Peptide_Data(tol, IMS_mass, LIMMA_pep);
            if ~isempty(Filtered_DS)
                Filtered_DS = [LIMMA_pep(1,:); Filtered_DS];
                %% Write output to excel
                [~,filename_without_ext,~] = fileparts(filename_with_path);
                out_file = strcat(dirname,'/Output/', filename_without_ext, out_file_end_part);
                xlswrite(out_file, Filtered_DS);
            else
                %             fprintf(fileID, 'FC_grt_0_Pvalue_DC: No matching IMS mass was found for %s file on %s\n', D(i).name, datetime('now'));
                %             fprintf('FC_grt_0_Pvalue_DC: No matching IMS mass was found for %s file on %s\n', D(i).name, datetime('now'));
                fn = strsplit(filename_with_path, '/');
                fprintf(Log_File_Id, '%s No matching IMS mass was found for %s file on %s\n', Log_File_Tag, fn{end}, datetime('now'));
                fprintf('%s No matching IMS mass was found for %s file on %s\n', Log_File_Tag, fn{end}, datetime('now'));
            end
        end        
    end
end