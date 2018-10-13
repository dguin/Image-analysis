names = dir('*.mat');
path = pwd;

%figure out temperature changes from mCherry calibration
jump = xlsread('jump.xlsx');
date1 = strsplit(pwd,'\');
date2 = strsplit(date1{end},'_');
date = str2double(date2{1,1});
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
        disp(strcat('save to file: ',num2str(date), '_mEGFP_',HSP_names{1,hsptype},'_temps'));
        save(strcat(num2str(date), '_mEGFP_',HSP_names{1,hsptype},'_temps'),'temps_jump63x', 'temp_change');
    elseif ~isempty(pgktype) && isempty(hsptype)
        disp(strcat('save to file: ',num2str(date), '_mEGFP_PGK', num2str(pgktype), '_mCh_temps'));
        save(strcat(num2str(date), '_mEGFP_PGK', num2str(pgktype), '_mCh_temps'),'temps_jump63x', 'temp_change');
    else
        Is_control_FPonly = input('Is this a mEGFP only control experiment? 1=yes, 0=no  ');
        if Is_control_FPonly == 1
            disp(strcat('save to file: ',num2str(date), '_mEGFP_temps'));
            save(strcat(num2str(date), '_mEGFP_temps'),'temps_jump63x', 'temp_change');
        else
            disp(strcat('save to file: ',num2str(date), '_mEGFP_mCh_temps'));
            save(strcat(num2str(date), '_mEGFP_mCh_temps'),'temps_jump63x', 'temp_change');
        end
    end
else
    disp(strcat('save to file: ',num2str(date), '_mEGFP_PGK', num2str(pgktype), '_',HSP_names{1,hsptype},'_temps'));
    save(strcat(num2str(date), '_mEGFP_PGK', num2str(pgktype), '_',HSP_names{1,hsptype},'_temps'),'temps_jump63x', 'temp_change');
end


for j = 1:size(names,1)
    i1 = strfind(names(j).name, 'red');
    i2 = strfind(names(j).name, 'temps');
    i3 = strfind(names(j).name, 'cells');
    i4 = strfind(names(j).name, 'finallist');
    i5 = strfind(names(j).name, 'compiled');
    i6 = strfind(names(j).name, 'matlab');
    split_name = strsplit(names(j).name,'.');
    if isempty(i1) && isempty(i2) && isempty(i3) && isempty(i4) && isempty(i5) && isempty(i6)
        char(split_name(1))        
        cellfind_multijumps_withsegmentation_redchannel_hist_thresh_sfr(pwd, char(split_name(1)), 60);
    end
end