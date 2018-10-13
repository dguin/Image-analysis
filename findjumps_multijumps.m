function findjumps_multijumps( datapath, dataname )
%Opens partially analysed matlab workspaces and implements jump separation
%   MATLAB workspaces with cells first analyzed by
%   cellfind_multijumps_withsegmentation. leftTrace and rightTrace are then
%   separated into jumps. Jumps are separated by manual or automatic peak
%   picking.

%Open matlab workspace 
    cd(datapath); 
    load(dataname, 'rightTrace', 'leftTrace', 'daTrace', 'time');
    split_name = strsplit(dataname,'.');
    split_name = char(split_name(1));
    %intitalize jump_traces_green to check sizes
    jump_traces_green = {};
    %Intialize variables
    numObjects = size(leftTrace,1);
    numJumps = 10;
    included_cells = zeros(numObjects, 1);
    pks = zeros(numObjects, numJumps+2);
    pks(:,1) = 16;
    pks(:,end) = 4900;
    for i = 1:numObjects
        fig_name = strsplit(dataname, '_cells');
        f1=openfig(strcat(char(fig_name(1)),'_aligned'));
        f = figure; 
        %plot differential of mean green channel signal. Points at which the
        %temperature is jumped should appear as spikes.
        plotyy(1:1:5759,diff(rightTrace(i,1:end)),1:1:5760,rightTrace(i,1:end)); 
        axis([0 inf -inf 1]);
        %skip cell is jumps don't look clean or aren't apparent as spikes.
        skip_loop = input('Enter 0 to continue, 1 to skip current cell ');
        close(f1);
        close(f);
        if skip_loop == 1
            continue
        else
            included_cells(i,1) = i;
        end
        %decide if peaks are picked manually or automatically
        pickhow = input('Pick peak manually? Y/N ', 's');
        %if no user input choose manual by default
        if isempty(pickhow)
            pickhow = 'N';
        end

        %for manual picking open figure of 1st derivative of green channel and
        %then let user pick out 10 points
        if (pickhow == 'Y' || pickhow == 'y')
            disp('Pick peaks now.');
            f = figure;
            plotyy(1:1:5759,diff(rightTrace(i,1:end)),1:1:5760,rightTrace(i,1:end));
            axis([0 inf -inf 1]);
            [pks_picked, ~] = getpts(f);
            %If X is a vector of length m, then Y = diff(X) returns a vector of
            %length m-1. Add 1 to the point indices to offset the index by 1 to
            %account for the loss.
            pks_picked = pks_picked + 1;
            %Make sure that the index picked is actually at the peak of the
            %spike.
            for k = 1:size(pks_picked,1)
                %convert the index to whole number and expand search for maxima
                %to 10 points centered around the peak picked
                peak_index = floor(pks_picked(k));
                [~,pks_picked(k)] = max(-diff(rightTrace(i,peak_index-30:peak_index+30)));
                pks_picked(k) = pks_picked(k) + peak_index - 30;
            end
            pks(i,2:size(pks_picked)+1) = pks_picked;
            close(f);

            faketrace = -diff(rightTrace(i,:));
            x=1:length(faketrace);
            x(pks_picked-1)=[];
            faketrace(x)=0;
            f = figure;
            subplot(1,2,1)
            findpeaks(faketrace,'Npeaks',numJumps,'MinPeakDistance', 300);
            hold on
            plot(-diff(rightTrace(i,:)),'k--');
            subplot(1,2,2)
            plot(rightTrace(i,:), 'k')
            hold on
            plot(pks_picked, rightTrace(i,pks_picked), 'ro');
            disp('If jumps selection looks good click on figure to proceed');
            waitforbuttonpress;
            close(f);
            clear faketrace x
        else
            f = figure;
            plot(diff(rightTrace(i,:)));
            axis([0 inf -inf 1]);
            grid on;
            %automatic recognition picks peaks according to a minimum height
            %cutoff. Only use if all spikes are over a minimum height cutoff
            %and are separated by atleast 450 points.
            cutoff = input('Enter cutoff height for automatic jump recognition: ');
            [~,pks_picked] = findpeaks(-diff(rightTrace(i,:)),'Npeaks',numJumps,'MinPeakDistance', 400, 'MinPeakHeight', cutoff);
            pks_picked = pks_picked+1;
            close(f);

            f = figure;
            subplot(1,2,1)
            findpeaks(-diff(rightTrace(i,:)),'Npeaks',numJumps,'MinPeakDistance', 400, 'MinPeakHeight', cutoff);
            hold on
            plot(-diff(rightTrace(i,:)), 'k--');
            subplot(1,2,2)
            plot(rightTrace(i,:), 'k')
            hold on
            plot(pks_picked, rightTrace(i,pks_picked), 'ro');
            disp('If jumps selection looks good click on figure to proceed');
            waitforbuttonpress;
            close(f);
            pks(i,2:size(pks_picked,2)+1) = pks_picked;
        end
        pks(i,:) = sort(pks(i,:),2);
        for j = 1:size(pks,2)-1 %separate video into jumps
            jump_traces_green{j,i} = rightTrace(i, pks(i,j) - 4: pks(i,j+1) - 5);
            jump_traces_red{j,i} = leftTrace(i, pks(i,j) - 4: pks(i,j+1) - 5);
            jump_traces_da{j,i} = daTrace(i, pks(i,j) - 4: pks(i,j+1) - 5);
            jump_time{j} = time(1, pks(i,j) - 4: pks(i,j+1) - 5);
        end
    clear pks_picked
    end

    if size(jump_traces_green,1)~=0           
        disp(strcat('save and append to file: ',split_name));
        save(split_name, 'jump_traces_green', 'jump_traces_red', 'jump_traces_da', 'jump_time', 'pks', 'included_cells', '-append');
    else
        save(split_name, 'included_cells', '-append');
    end

end

