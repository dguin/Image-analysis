function Compile_all_cells(folders)
    % intialize following variables
    % folders{1:n,1} = foldernames;
    % Then run --> Get_cell_DA_and_Red;
    % for ii = 1:size(folders,1)
    allcells_accepted = [];
    allcells_accepted_fretE = [];
    rTrace = [];
    dTrace = [];
    slopeind_end = 5;
    slopeind_start = 1;
    current_path = 'C:\Users\dguin\Documents\MATLAB\Hsp70\In-cell data\20180717_incelldata_reanalysis';
    for ii = 1:size(folders,1)
        cd(folders{ii,1});
        disp(folders{ii,1});
        date1 = strsplit(pwd,'\');
        date2 = strsplit(date1{end},'_');
        date = str2double(date2{1,1});
        jump = xlsread('jump.xlsx');
        jump63x;
        HSP_names = cell(1,4);
        HSP_names{1,1} = 'mHSP';
        HSP_names{1,2} = 'mHSPK71M';
        HSP_names{1,3} = 'mHSC';
        HSP_names{1,4} = 'mHSCK71M';
        Is_control = input('Is this a control experiment? 1=yes, 0=no  ');
        pgktype = input('Enter PGK variant number if available, 0-3, press enter to leave empty  ');
        hsptype = input('Enter HSP variant(1-mHSP70, 2-mHSP70K71M, 3-mHSC70, 4-mHSC70K71M), press enter to leave empty  ');
        if Is_control == 1
            if isempty(pgktype) && ~isempty(hsptype)
                string_dataname = strcat(num2str(date), '_mEGFP_',HSP_names{1,hsptype},'_compiled');
                string_allcellscompiled_name = strcat('mEGFP_',HSP_names{1,hsptype},'_compiled');
            elseif ~isempty(pgktype) && isempty(hsptype)
                string_dataname = strcat(num2str(date), '_mEGFP_PGK', num2str(pgktype), '_mCherry_compiled');
                string_allcellscompiled_name = strcat('mEGFP_PGK', num2str(pgktype), '_mCherry_compiled');
            else
                Is_control_FPonly = input('Is this a mEGFP only control experiment? 1=yes, 0=no  ');
                if Is_control_FPonly == 1                    
                    string_dataname = strcat(num2str(date), '_mEGFP_compiled');
                    string_allcellscompiled_name = strcat('mEGFP_compiled');
                else
                    string_dataname = strcat(num2str(date), '_mEGFP_mCherry_compiled');
                    string_allcellscompiled_name = strcat('mEGFP_mCherry_compiled');
                end
            end
        else
            string_dataname = strcat(num2str(date), '_mEGFP_PGK', num2str(pgktype),'_',HSP_names{1,hsptype},'_compiled');
            string_allcellscompiled_name = strcat('mEGFP_PGK', num2str(pgktype),'_',HSP_names{1,hsptype},'_compiled');
        end

        % dataname = strcat(string_dataname, num2str(pgktype), '_', num2str(date));
        % dataname = strcat(string_dataname, '_', num2str(date));
        includeinmean = [];
        erasedfrommean = [];
        da_interp = [];
        da_interp_norm = [];
        l = 1;
        load('finallist.mat')
        for j = 1:size(finalcelllist,1)
            load(finalcelllist(j).name, 'jump_traces_da', 'jump_traces_green', 'jump_traces_red', 'jump_traces_area', 'red_mean', 'rightTrace', 'leftTrace', 'daTrace');
            disp(finalcelllist(j).name);
            if isempty(jump_traces_da{1,finalcelllist(j).index})
                disp('passed')
                continue
            end
            dataname(j,:) = strsplit(finalcelllist(j).name, '.');
            %dataname{j,2} = [];
            for i = 1:11
                xda = jump_traces_da{i,finalcelllist(j).index};
                xg = jump_traces_green{i,finalcelllist(j).index};
                xr = jump_traces_red{i,finalcelllist(j).index};
                xa = jump_traces_area{i,finalcelllist(j).index};
                if (i == 1 || i == 11)  
                    da(i,l) = mean(xda(10:end-10));
                    green(i,l) = mean(xg(10:end-10));
                    red(i,l) = mean(xr(10:end-10));
                    cellarea(i,l) = mean(xa(10:end-10));
                else
                    da(i,l) = mean(xda(end-120:end-20));
                    green(i,l) = mean(xg(end-120:end-20));
                    red(i,l) = mean(xr(end-120:end-20));
                    cellarea(i,l) = mean(xa(end-120:end-20));
                end
            end
            for i = 1:10
                effective_mEGFP_temp_change(i,l) = ((green(i,l) - green(i+1,l))./green(i,l))*(100/1.63);
            end
            average_red = mean(red_mean(finalcelllist(j).index,:),2);
            amber_Cell_Red(l,1) = average_red;
            rTrace(l,:) = rightTrace(finalcelllist(j).index,:);
            dTrace(l,:) = daTrace(finalcelllist(j).index,:);
            lTrace(l,:) = leftTrace(finalcelllist(j).index,:);
            l = l + 1;
        end
        dataname(:,2) = [];
        temperatures = temps_jump63x;
        disp(temperatures(end));
        Final_temp = input('Enter final temperture for interpolation, must be even ');
        temps_interp = (20:2:Final_temp)';
        filter_Cells_Red_Green; 
        disp(strcat(num2str(size(filtered_da,2)), 'cells passed'));
        for i = 1:size(filtered_da,2)
%             da_interp(:,i) = pchip(temperatures, smooth(filtered_da(:,i)), temps_interp);
%             da_interp_norm(:,i) = da_interp(:,i)./da_interp(1,i);
%             slopes = polyfit(temps_interp(1:4), da_interp(1:4,i), 1);
%             slope_corrected_da = da_interp(:,i) - slopes(1).*(temps_interp - temps_interp(1));
%             slopes_norm = polyfit(temps_interp(1:4), da_interp_norm(1:4,i), 1);
%             slope_corrected_da_norm = da_interp_norm(:,i) - slopes_norm(1).*(temps_interp - temps_interp(1));
            
            fig_name = strsplit(dataname{i,1}, '_cells');
            f1=openfig(strcat(char(fig_name(1)),'_aligned'));

            cellarea_interp(:,i) = pchip(temperatures, (filtered_cellarea(:,i)), temps_interp);
            green_interp(:,i) = pchip(temperatures, (filtered_green(:,i)), temps_interp);
            green_interp_norm(:,i) = green_interp(:,i)./green_interp(1,i);
            slopes_green = polyfit(temps_interp(slopeind_start:slopeind_end), green_interp(slopeind_start:slopeind_end,i), 1);
            slope_corrected_green = green_interp(:,i) - slopes_green(1).*(temps_interp - temps_interp(1));
            slopes_norm_green = polyfit(temps_interp(slopeind_start:slopeind_end), green_interp_norm(slopeind_start:slopeind_end,i), 1);
            slope_corrected_norm_green = green_interp_norm(:,i) - slopes_norm_green(1).*(temps_interp - temps_interp(1));
            
            red_interp(:,i) = pchip(temperatures, (filtered_red(:,i)), temps_interp);
            red_interp_norm(:,i) = red_interp(:,i)./red_interp(1,i);
            slopes_red = polyfit(temps_interp(slopeind_start:slopeind_end), red_interp(slopeind_start:slopeind_end,i), 1);
            slope_corrected_red = red_interp(:,i) - slopes_red(1).*(temps_interp - temps_interp(1));
            slopes_norm_red = polyfit(temps_interp(slopeind_start:slopeind_end), red_interp_norm(slopeind_start:slopeind_end,i), 1);
            slope_corrected_norm_red = red_interp_norm(:,i) - slopes_norm_red(1).*(temps_interp - temps_interp(1));
            
            da_interp = green_interp./red_interp;
            fretE_interp = red_interp./(red_interp+green_interp);
            slope_da_interp = polyfit(temps_interp(slopeind_start:slopeind_end), da_interp(slopeind_start:slopeind_end,i), 1);
            slope_corrected_da_interp = da_interp(:,i) - slope_da_interp(1).*(temps_interp - temps_interp(1));
            slope_corrected_da_interp_norm(:,i) = slope_corrected_da_interp./slope_corrected_da_interp(1);
            da_interp_norm(:,i) = da_interp(:,i)./da_interp(1,i);
            fretE_interp_norm(:,i) = fretE_interp(:,i)./fretE_interp(1,i);
            slope_corrected_da = slope_corrected_green./slope_corrected_red;

            f = figure;
            title(strcat('Cell', num2str(i)));
            subplot(2,3,1)       
            plot(temperatures, filtered_da(:,i),'bo');
            hold on
            plot(temps_interp, da_interp(:,i),'b-');
            plot(temps_interp, slope_corrected_da_interp,'kd-');
            lgd1 = legend('exp data', 'interpolated data', 'slope corrected interp');
            lgd1.FontSize = 4;
            yL = get(gca,'YLim');
            line([46 46],yL,'Color','r');
            subplot(2,3,2)
            plot(temps_interp,slope_corrected_da./slope_corrected_da(1),'rs-');
            hold on
            plot(temps_interp,slope_corrected_da_interp_norm(:,i), 'm^--');
            plot(temps_interp,1./(1+(slope_corrected_da_interp_norm(:,i)-1)), 'bX--');
            lgd2 = legend('corrected red-green', 'corrected da');
            lgd2.FontSize = 4;
            axis([20 50 -inf inf]);
            yL = get(gca,'YLim');
            line([46 46],yL,'Color','r');
            subplot(2,3,3)
            plot(1:1:5759, diff(rTrace(i,:)));
            axis([0 5800 -inf inf]);
            subplot(2,3,4)
            plot(diff(dTrace(i,:)));
            axis([0 5800 -inf inf]);
            subplot(2,3,5)
            plot(1:1:5760, rTrace(i,:));
            axis([0 5800 -inf inf]);
            subplot(2,3,6)
            plot(temps_interp, cellarea_interp(:,i))
            axis([20 50 -inf inf]);
            disp(strcat('Saving traces to file: ', char(dataname{i,1}), '_', num2str(i), '_traces.jpg'));
            skip_loop = input('Enter 0 to accept cell to calculate mean, 1 to skip current cell ');
            saveas(f, strcat(char(dataname{i,1}), '_', num2str(i),'_traces.jpg'));
            if skip_loop == 1
                erasedfrommean = cat(1, erasedfrommean, i); 
                close(f);
                close(f1)
                cells_kept_greenred_threshholding(1,i) = 0;
            else
                includeinmean = cat(1, includeinmean, i);
                close(f);
                close(f1)
            end
        end
        allcells_accepted = horzcat(allcells_accepted,da_interp_norm(:,includeinmean));
        allcells_accepted_fretE = horzcat(allcells_accepted_fretE,fretE_interp_norm(:,includeinmean));
        %allcells_accepted = horzcat(allcells_accepted,da_interp_norm(:,includeinmean));
        cells_kept_greenred_threshholding(:,sum(cells_kept_greenred_threshholding,1)==0) = []; %remove empty columns
        cd(current_path);
        save(string_dataname,'da','red','green','amber_Cell_Red','filtered_da', 'filtered_green', 'filtered_red', 'removed_da', 'removed_green', 'removed_red', 'temperatures', 'date', 'da_interp', 'temps_interp', 'da_interp_norm', 'includeinmean', 'erasedfrommean','Green_Red','cells_kept_greenred_threshholding', 'effective_mEGFP_temp_change','amber_cell_red_filtered', 'dataname','slope_corrected_da_interp_norm','slopeind_end', 'slopeind_start');
        clear da red green amber_Cell_Red filtered_da filtered_green filtered_red removed_da removed_green removed_red temperatures date da_interp temps_interp da_interp_norm includeinmean erasedfrommean Green_Red cells_kept_greenred_threshholding rTrace dTrace dataname
    end
    temps_interp = 20:2:Final_temp;
    assignin('base', 'temps_interp', temps_interp);
    mean_da_interp = mean(allcells_accepted, 2);
    assignin('base', 'mean_da_interp', mean_da_interp);
    std_da_interp = std(allcells_accepted, 1, 2);
    assignin('base', 'std_da_interp', std_da_interp);
    stdmean_da_interp = std_da_interp/sqrt(size(allcells_accepted,2));
    assignin('base', 'stdmean_da_interp', stdmean_da_interp);
    slope_mean_da_interp = polyfit(temps_interp(slopeind_start:slopeind_end)', mean_da_interp(slopeind_start:slopeind_end),1);
    assignin('base', 'slope_mean_da_interp', slope_mean_da_interp);
    slope_corrected_mean_DA_interp = mean_da_interp - slope_mean_da_interp(1).*(temps_interp - temps_interp(1))';
    assignin('base', 'slope_corrected_mean_DA_interp', slope_corrected_mean_DA_interp);
    
    mean_fretE_interp = mean(allcells_accepted_fretE, 2);
    assignin('base', 'mean_fretE_interp', mean_fretE_interp);
    std_fretE_interp = std(allcells_accepted_fretE, 1, 2);
    assignin('base', 'std_da_interp', std_fretE_interp);
    stdmean_fretE_interp = std_fretE_interp/sqrt(size(allcells_accepted_fretE,2));
    assignin('base', 'stdmean_da_interp', stdmean_fretE_interp);
    slope_mean_fretE_interp = polyfit(temps_interp(slopeind_start:slopeind_end)', mean_fretE_interp(slopeind_start:slopeind_end),1);
    assignin('base', 'slope_mean_da_interp', slope_mean_fretE_interp);
    slope_corrected_mean_fretE_interp = mean_fretE_interp - slope_mean_fretE_interp(1).*(temps_interp - temps_interp(1))';
    assignin('base', 'slope_corrected_mean_DA_interp', slope_corrected_mean_fretE_interp);
    
    save(string_allcellscompiled_name,'allcells_accepted','temps_interp','mean_da_interp','std_da_interp','stdmean_da_interp','slope_mean_da_interp','slope_corrected_mean_DA_interp','allcells_accepted_fretE','mean_fretE_interp','std_fretE_interp','stdmean_fretE_interp','slope_mean_fretE_interp','slope_corrected_mean_fretE_interp','-v7.3')
end