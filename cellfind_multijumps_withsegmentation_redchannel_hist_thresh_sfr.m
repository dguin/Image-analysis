function cellfind_multijumps_withsegmentation_redchannel_hist_thresh_sfr( datapath, dataname, fps, thresh )
    %cellfind_multijumps_withsegmentation identifies cells using the
    %segmentation protocol in 
    %https://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
    %   1. Threshold cells and align channels
    %   2. Save partially processed data
    %   3. Convert labels to .avi movie for easier visualization of area
    %   changes

    %load data 
    strpath = pwd; 
    cd(datapath); 
    load(dataname, '-mat'); 
    
    red_filename = strcat(dataname, 'red.mat');
    red_exist = exist(red_filename);
    
    mincellarea = 1000;
    constant_thresh = thresh;
    
    %Get video info 
    [height, width, numFrames] = size(dataright); 
    time = (1:numFrames)*1/fps - 1/fps; 

    %first threshold green channel(dataright) using frame 1 and pick out cell
    %red channel image are low contrast so use adapthisteq to increase 
    %contrast slightly. adapthisteq implements a technique called 
    %contrast-limited adaptive histogram equalization, or CLAHE. 
    %im2bw converts the image to logical threshholded image
    %thresholded image is then fill, dilated and cleared of amll objects 
    %<100 pixels to avoid getting small cells or areas of noise that make 
    %it throught the thresholding algorithm.
    alignfr = 1;
    I_eq = adapthisteq(uint8(dataleft(:,:,alignfr)), 'Distribution', 'rayleigh');
    
    %Before proceeding align both channels
    %Estimate the transformation needed to align the images using 
    %imregtform. Imregtform estimates geometric transformation that aligns 
    %two 2-D or 3-D images
    [optimizer, metric] = imregconfig('multimodal');
    tformEstimate = imregtform(uint8(dataright(:,:,1)), I_eq, 'rigid', optimizer, metric);
    
    %Apply estimated geometric transform to the moving image. This example 
    %uses the 'OutputView' parameter to obtain a registered image the same 
    %size and with the same world limits as the reference image.
    dataright_realigned = imwarp(uint8(dataright),tformEstimate,'OutputView',imref2d(size(I_eq)));

    f1=figure;
    %images overlayed
    Im1 = imfuse(I_eq, dataright_realigned(:,:,1), 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    %images side-by-side
    Im2 = imfuse(I_eq, dataright_realigned(:,:,1), 'montage');
    %title('Panel 1 images overlayed');
    
    %initial threshold frame
    threshfr = 1;
    
    %panel 2 graythresh with non adaptive processing
    subplot(3,3,2);
    for i=1:100
        setThresh(i)=(i-1)*(1/100);
        getThresh_green(i)=sum(sum(im2bw(uint8(dataright_realigned(:,:,threshfr)),setThresh(i))));
        getThresh_red(i)=sum(sum(im2bw(uint8(dataleft(:,:,threshfr)),setThresh(i))));
    end
    plotyy((((1:100)-1)*(1/100)), getThresh_green, (((1:100)-1)*(1/100)), getThresh_red);
    hold on
    plot(graythresh(dataright_realigned(:,:,threshfr))*ones(100,1),getThresh_green,'b');
    plot(graythresh(uint8(dataleft(:,:,threshfr)))*ones(100,1),getThresh_red,'r');
    xlabel('Thresholds');
    ylabel('Pixels included');
    %title('Panel 2 graythresh with non adaptive processing');
    
    [labels, cellStats, thresh, numObjects, NP] = thresh_frames(uint8(dataleft(:,:,threshfr)), constant_thresh, mincellarea, numFrames);
    
    %panel 3 graythresh mask with non adaptive processing
    subplot(3,3,3);
    BW = (labels(:,:,1).*dataright_realigned(:,:,threshfr))*255;
    %[~,temp_image] = image_binarize(BW, mincellarea);
    imshowpair(BW, dataright_realigned(:,:,threshfr), 'montage');
    %title('Panel 3 graythresh cell mask with non adaptive processing');
    clear temp_image BW
    
    for i = 1:numObjects
        %Convert label matrix into binary for calculating cell size
        %convert labels to uint8 and write values in labels_cellsize
        labels_cellsize(:,:,i) = uint8(labels);
        %convert all non zero values in labels_movie to 1 to calculate cell
        %size
        labels_cellsize(labels~=0 & labels~=i) = 1;
    end
    
    clear i
    
    %panel 5 threshold values vs. frame
    subplot(3,3,5)
    bar(1:1:numObjects, thresh, 'k');
    hold on
    bar(numObjects+1, NP, 'r');
    xlim([0 numObjects+2])
    xlabel('Cell Number');
    ylabel('Threshold value');
    
    %panel 6 bar distribution of a red channel image
    subplot(3,3,6)
    [counts,edges] = histcounts(reshape(dataleft(:,:,1),1,size(dataleft,1)*size(dataleft,2)),256,'binmethod','integers');
    centers = (edges(1:end-1) + edges(2:end))/2;
    plot(fit(centers', counts','gauss2'), centers, counts, 'b')
    legend('off');
    xlabel('Pixel Intensity');
    ylabel('Number of pixels');
    
    %initialize variables
    leftTrace = zeros(numObjects, numFrames);
    rightTrace = zeros(numObjects, numFrames);
    daTrace = zeros(numObjects, numFrames);
    cell_size = zeros(numFrames, numObjects);
    
    for j = 1:numObjects
        %find row/column index assigned to cell 'j'
        [currentRow, currentColumn] = find(labels==j);
        if isempty(currentRow)
            continue
        end
        for i = 1:numFrames
            %pick out intensity values for pixels in currentRow and
            %currentColumn
            for k = 1:length(currentRow)
                currentLeft(k,1) = dataleft(currentRow(k), currentColumn(k), i);
                currentRight(k,1) = dataright_realigned(currentRow(k), currentColumn(k), i);                
            end  
            if isempty(currentLeft)
                continue
            end
            [low_left] = find(abs(currentLeft)<NP);
            currentLeft(low_left) = [];
            currentRight(low_left) = [];
            cell_size(i,j)=size(currentLeft,1);
            %Average cell pixels leftChannel and rightChannel
            leftTrace(j,i) = mean(double(currentLeft));
            rightTrace(j,i) = mean(double(currentRight));
            daTrace(j,i) = rightTrace(j,i)./leftTrace(j,i);
            %save currentLeft and currentRight in cells Left and Right
            %Left{i,j} = currentLeft;
            %Right{i,j} = currentRight;
            clear currentLeft currentRight
        end
    end
    
    if ~isempty(find(leftTrace==0,1))
        leftTrace(leftTrace == 0) = [];
        rightTrace(rightTrace == 0) = [];
        daTrace(daTrace == 0) = [];

        leftTrace = reshape(leftTrace, size(leftTrace,2)/5760, 5760);
        rightTrace = reshape(rightTrace, size(rightTrace,2)/5760, 5760);
        daTrace = reshape(daTrace, size(daTrace,2)/5760, 5760);
    end
    
    %panel 4 cell size vs. frame
    subplot(3,3,4)
    plot(1:1:numFrames, cell_size);
    xlabel('Cell Number');
    ylabel('Pixels included');
    xlim([0 6000])
    
    clear currentRow currentColumn i j k
    
    %panel 6 leftTrace values vs. time
    subplot(3,3,7)
    plot(time, leftTrace, 'r');
    xlabel('time (secs)');
    ylabel('Red Trace');
    
    %panel 7 leftTrace values vs. time
    subplot(3,3,8)
    plot(time, rightTrace, 'g');
    xlabel('time (secs)');
    ylabel('Green Trace');
    
    %panel 8 daTrace values vs. time
    subplot(3,3,9)
    plot(time, daTrace, 'b');
    xlabel('time (secs)');
    ylabel('D/A Trace');
    
    %VideoName will consist of slide and cell number and is saved in parent
    %datapath
    VideoName = strsplit(dataname, '_');
    %example Video datapath ---> datapath\s1c1b.avi
    LabelVideo = VideoWriter(fullfile(datapath, strcat(char(VideoName{1,end}),'.avi')));
    LabelVideo.FrameRate = 60;
    open(LabelVideo)
    disp(strcat('Writing label file to video ->', strcat(char(VideoName{1,end}),'.avi')));
    for i = 1:numFrames
        %Convert label matrix into movie for manual checking
        %convert labels to uint8 and write values in labels_movie
        labels_movie_green = labels.*dataright_realigned(:,:,i);
        labels_movie_red = labels.*uint8(dataleft(:,:,i));
        %Increase contrast of images
        labels_contrast_inc_green = (labels_movie_green);
        labels_contrast_inc_red = (labels_movie_red);
        %Fuse both constrast adjusted images together
        labels_movie = imfuse(labels_contrast_inc_green, labels_contrast_inc_red, 'montage');
        %write labels_movie into video
        writeVideo(LabelVideo,labels_movie);
        clear labels_movie_green labels_movie_red labels_contrast_inc_green labels_contrast_inc_red
    end
    close(LabelVideo)
    
    disp(strcat('save all variables to file:  ',dataname,'_cells.mat'));
    save(strcat(dataname, '_cells'),'labels', 'cellStats', 'time', 'leftTrace', 'rightTrace', 'daTrace', 'height', 'width', 'thresh', '-v7.3');
    
    %Handle amber only red channel
    if red_exist==2
        load(red_filename);
        %Align amber images
        tformEstimate_red = imregtform(uint8(red(:,:,1)), I_eq, 'rigid', optimizer, metric);
        red_realigned = imwarp(uint8(red),tformEstimate_red,'OutputView',imref2d(size(I_eq)));
        %panel 3 Amber red overlayed with FRET red
        Imred = imfuse(I_eq, red_realigned(:,:,1), 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
        %display overlayed FRET red with green and overlayed with amber red
        %side by side
        Im3 = imfuse(Im1, Imred, 'montage');
        %display all images side by side
        Im4 = imfuse(Im3, Im2, 'montage');
        subplot(3,3,1);
        imshow(Im4);
        %save aligned images
        saveas(f1, strcat(dataname, '_aligned.jpg'));
        saveas(f1, strcat(dataname, '_aligned.fig'));
        close(f1);
        for j = 1:numObjects
            [currentRow_red, currentColumn_red] = find(labels(:,:,1)==j);
            if isempty(currentRow_red)
                continue
            end            
            for k = 1:length(currentRow_red)
                currentRed(k,:,j) = red_realigned(currentRow_red(k), currentColumn_red(k),:);                        
            end
            red_mean(j,:) = mean(double(currentRed(:,:,j)));
        end
        if ~isempty(find(red_mean==0,2))
            red_mean(red_mean == 0) = [];    
            red_mean = reshape(red_mean, size(red_mean,2)/60, 60);
        end
        disp(strcat('save amber red to file: ',dataname,'_cells.mat'));
        save(strcat(dataname, '_cells'), 'red_mean', '-append');
    else
        %display overlayed FRET red with green and side by side
        Im3 = imfuse(Im1, Im2, 'montage');
        subplot(3,3,1);
        imshow(Im3);
        %save aligned images
        disp(strcat('No red video for ->', dataname));
        saveas(f1, strcat(dataname, '_aligned.jpg'));
        saveas(f1, strcat(dataname, '_aligned.fig'));
        close(f1);
    end
    clear I_eq
end

function [labels, cellStats, thresh, numObjects, NP] = thresh_frames(data, constant_thresh, mincellarea, ~)
    %temporarily convert images to binary using graythresh
    BW0 = imadjust(data(:,:,1));
    BW = im2bw(BW0, graythresh(BW0)); 
    clear BW0
    [cells,~] = image_binarize(BW, mincellarea);
    cell_IDs = regionprops(cells, 'BoundingBox');
    cell_Labels = labelmatrix(cells);
    numObjects = size(cell_IDs, 1);
%     numObjects = 1;
    frames_to_thresh = 1;
    thresh = zeros(numObjects,frames_to_thresh);
    labels = uint8(zeros(size(data,1),size(data,2),frames_to_thresh));       
                
    for j = 1:numObjects
        I_cell = data;
        I_rect = [1 1 size(data,2)-1 size(data,1)-1];
        width_cropim = round(I_rect(1,1)+I_rect(1,3));
        height_cropim = round(I_rect(1,2)+I_rect(1,4));        
        for i = 1:frames_to_thresh
            %First fit the image using single gaussian
            %The histogram of the values from 0 to 255 is stored in counts,
            %Edges gives the left and right edge of the bin
            %For instance, count(1) contains the number of occurrence of the 
            %value zero in the image.
            [counts,edges] = histcounts(reshape(I_cell(:,:,i),1,size(I_cell,1)*size(I_cell,2)),256,'binmethod','integers');
            %centers calculates the mid point of each bin; should be te integer
            %pixel intensity value between 0-255
            centers = (edges(1:end-1) + edges(2:end))/2;
            %fit the histogram using a Gauss1 fit. This should mainly fit the
            %large peak due to background pixels. 
            Gauss_noise = fit(centers', counts', 'gauss1'); 
            %Find out the noise pixel value mu+sigma
            NP = round(Gauss_noise.b1 + (Gauss_noise.c1/sqrt(2)));

            if numObjects>1 && i==1
                clear I_cell I_rect
                %Choose smallest box containing foreground pixels
                I_rect = cell_IDs(j).BoundingBox;
                I_rect(1,1:2) = ceil(I_rect(1,1:2)-10);
                I_rect(1,3:4) = I_rect(1,3:4)+10;     
                width_cropim = round(I_rect(1,1)+I_rect(1,3));
                height_cropim = round(I_rect(1,2)+I_rect(1,4));
                if ~isempty(find(I_rect<1))
                    I_rect(find(I_rect<1)) = 1;
                end
                if I_rect(1,2)+height_cropim>size(data,1)
                    height_cropim = size(data,1) - 1;
                end
                if I_rect(1,1)+width_cropim>size(data,2)
                    width_cropim=size(data,2) - 1;
                end
                %Make sure number of cells is one in I_cell
                cell_IDs_afterautocrop = cell_Labels(I_rect(1,2):height_cropim, I_rect(1,1):width_cropim);
                [~,~,numObjects_afterautocrop] = find(cell_IDs_afterautocrop~=0 & cell_IDs_afterautocrop~=j);
                if ~isempty(numObjects_afterautocrop)
                    %if still more than one cell in I_cell do manual
                    %cropping
                    clear row_cellj col_cellj cell_IDs_afterautocrop numObjects_afterautocrop
                    f2 = figure;
                    image_to_show = imadjust(data(:,:,1));
                    imshow(image_to_show);
                    I_rect_new = imrect(gca, I_rect);
                    I_rect = wait(I_rect_new);
                    close(f2);
                end
                I_cell = data(I_rect(1,2):height_cropim, I_rect(1,1):width_cropim, :);
            end
            %Fit the the data again using calculated fit parameters in
            %Gauss_noise
            fity = Gauss_noise(centers);
            %Find the first index where fity is less than 1 after maximum
            %of Gauss_noise has occured. Assume this where the second gaussian
            %with foreground signal distribution starts
            [~,max_fity_idx] = max(fity);
            row_ind = find(fity(max_fity_idx:end)<(0.01*Gauss_noise.a1),1) + max_fity_idx - 1;
            counts_interp = interp1(centers(row_ind:end), counts(row_ind:end), centers(row_ind):0.2:centers(end));
            %Fit the second smaller gaussian with brighter pixels starting at
            %row_ind
            if constant_thresh == 0    
                GaussB = fit((centers(row_ind):0.2:centers(end))', counts_interp', 'gauss1');
                %Find threshold at left side FWHM x value for GaussB
                %specifically at b1-(c1/sqrt(2)) or mu-sigma. Choose only top
                %25% values
                thresh(j,i) = floor(GaussB.b1 - (0.25*(GaussB.c1/sqrt(2))));
            else
                thresh(j,i) = constant_thresh;
            end
            if thresh(j,i) < centers(row_ind)
                Gauss_noise_interp = Gauss_noise(centers(1):0.2:centers(end));
                [~,max_Gauss_noise_interp_idx] = max(Gauss_noise_interp);
                row_ind_thresh_cutoff = find(Gauss_noise_interp(max_Gauss_noise_interp_idx:end)<(0.01*Gauss_noise.a1),1) + max_Gauss_noise_interp_idx - 1;
                disp('Threshold too low, using gaussian noise cutoff')
                thresh(j,i) = round(centers(1) + (row_ind_thresh_cutoff-1)*0.2);
            end
            disp(thresh(j,i));
            %start image binarization based on thresh(i)
            BW = I_cell(:,:,i);
            BW(BW<thresh(j,i))=0;
            BW = logical(BW);
            [cells,~] = image_binarize(BW, mincellarea);
            CellNum_ID(:,:,i) = labelmatrix(cells);
            if max(max(CellNum_ID(:,:,i)))>1
                temp_cellstats = regionprops(cells, 'Area');
                temp_cellSizes = cat(1, temp_cellstats.Area);
                [~,idx_bigcell] = max(temp_cellSizes);
                temp_CellNum_ID = CellNum_ID(:,:,i);
                temp_CellNum_ID(temp_CellNum_ID~=idx_bigcell) = 0;
                temp_CellNum_ID(temp_CellNum_ID==idx_bigcell) = 1;
                CellNum_ID(:,:,i) = temp_CellNum_ID;
            end
            clear counts edges centers Gauss_noise fity GaussB row_ind BW cells temp_cellSizes temp_CellNum_ID temp_cellstats
        end
        CellNum_ID(CellNum_ID~=0) = CellNum_ID(CellNum_ID~=0) + j - 1;
        labels(I_rect(1,2):height_cropim, I_rect(1,1):width_cropim, :) = CellNum_ID(:,:,:);
        clear I_cell I_rect CellNum_ID height_cropim width_cropim
    end  
    cellStats{:,1} = regionprops(labels(:,:,1), data, 'BoundingBox', 'Centroid', 'PixelList', 'MeanIntensity', 'Area');
    clear i j
end

function [cells,BW6] = image_binarize(BW, mincellarea)
    %BW = im2bw(data, thresh_im);
    BW1 = bwareaopen(BW, mincellarea);
    BW2 = bwperim(BW1);
    BW3 = imfill(BW2,'holes');
    BW4 = bwmorph(BW3, 'clean');
    seD = ones(3, 'uint8');
    BW5 = imclose(BW4, seD);
    BW6 = bwareaopen(BW5, mincellarea);
    %identify the cell using bwconncomp, bwconncomp recognizes connected
    %objects. Use bwconncomp component structure to create a matrix of 
    %cells where each cell is identified by a array of pixels labeled 
    %1,2,... Cell 1 pixels for example are labeled 1, cell 2 2 and so on
    cells =  bwconncomp(BW6); 
end

