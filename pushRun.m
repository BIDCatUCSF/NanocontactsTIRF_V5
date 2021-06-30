function pushRun(~,~,checkNewFolder,checkMean,checkMedian,checkAuto,checkManual,...
    checkBackground,checkGaussian,guiNC,handles)
%pushRun The 'RUN!' button essentially runs the program to analyze the
% nanocontacts. Includes most of the error handling. R2015b
% 
%
% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% May 2017 - June 2021
    %% Load the structure and assign variables
    tic
    structParameters = getappdata(guiNC,'structParameters');
    
    tcrName = structParameters.tcrName;
    maskName = structParameters.maskName;
    contourName = structParameters.contourName;
    newFolder = structParameters.newFolder;
    backgroundDisc = structParameters.backgroundDisc;
    gaussianSigma = structParameters.gaussianSigma;
    SDoverNoise = structParameters.SDoverNoise;
    
    timeRes = structParameters.timeRes;
    
    xMin = structParameters.xMin;
    xMax = structParameters.xMax;
    yMin = structParameters.yMin;
    yMax = structParameters.yMax;
    
    pathDir = strcat(structParameters.pathDir,'\');
    
    
    %% Some initial error handling
    % Check to make sure there was file input
    if strcmp(tcrName,'') == 1
        warndlg('No TCR file was selected','No File Found');
        return
    elseif strcmp(maskName,'') == 1
        warndlg('No Mask file was selected','No File Found');
        return
    elseif strcmp(contourName,'') == 1
        warndlg('No Contour file was selected','No File Found');
        return
    % Check to make sure an output directory exists
    elseif exist(pathDir,'dir') ~= 7
        warndlg('Output path does not exist','No Folder Found');
        return
    % Check if new folder already exists. If it does, offer the user to
    % overwrite the files, or choose a new folder name. If it does not,
    % create the new folder in the output path directory.
    elseif get(checkNewFolder,'Value') == 1
        if exist(strcat(pathDir,newFolder),'dir') == 7
            choice = questdlg('Warning','Analysis Directory','Overwrite','Cancel','Cancel');
            switch choice
                case 'Overwrite'
                    outFolder = strcat(pathDir,newFolder,'\');
                case 'Cancel'
                    return
            end % switch
        else
            mkdir(pathDir,newFolder);
            outFolder = strcat(pathDir,newFolder,'\');
        end % if
    else
    outFolder = pathDir;
    end % if
    
    % Make sure only one Signal Determination method is chosen
    if get(checkMean,'Value') == 1 && get(checkMedian,'Value') == 1
        warndlg('Please select only one Signal Determination Method','Select only one');
        return
    elseif get(checkMean,'Value') == 0 && get(checkMedian,'Value') == 0
        warndlg('Please select a Signal Determination Medthod','Select one');
        return
    end
    
    % Make sure only one ROI Selection method is chosen
    if get(checkAuto,'Value') == 1 && get(checkManual,'Value') == 1
        warndlg('Please select only one ROI Selection Method','Select only one');
        return
    elseif get(checkAuto,'Value') == 0 && get(checkManual,'Value') == 0
        warndlg('Please select one ROI Selection Method','Select one');
        return
    end
    
    %%
    structParameters.status = 'Status: Opening Files...';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
    % Open the files into stacks
    TCRStack = openTIFF(pathDir,tcrName);
    maskStack = openTIFF(pathDir,maskName);
    contourStack = openTIFF(pathDir,contourName);
    % Ensure the stacks are ones and zeros
    maskStack = makeBinary(maskStack,1);
    contourStack = makeBinary(contourStack,1);
    
    dimensions = size(TCRStack);
    % Invert maskStack
    maskStack2 = ones(dimensions);
    for ii = 1:dimensions(3)
        for jj = 1:dimensions(1)
            for kk = 1:dimensions(2)
                if maskStack(jj,kk,ii) == 1
                    maskStack2(jj,kk,ii) = 0;
                end % if
            end % if
        end % if
    end % if
    % Save the watershed stack to tif
    structParameters.status = 'Status: Creating Watershed...';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
    watershedMaskStack = watershedMask(maskStack2);
    outName = strcat(outFolder,'01_Watershed_',tcrName);
    writeTiff(watershedMaskStack,outName);
    
    %% User Inputs
    % Background subtract if asked for
    if get(checkBackground,'Value') == 1
        structuredArea = strel('disk',backgroundDisc);
        backgroundStack = zeros(dimensions);
        for ii = 1:dimensions(3)
            backgroundStack(:,:,ii) = imtophat(TCRStack(:,:,ii),structuredArea);
        end % for
    else
        backgroundStack = TCRStack;
    end % if
    outName = strcat(outFolder,'02_BackgroundSubtracted_',tcrName);
    writeTiff(backgroundStack,outName);
    % Gaussian filter is asked for
    if get(checkGaussian,'Value') == 1
        gaussianStack = zeros(dimensions);
        for ii = 1:dimensions(3)
            gaussianStack(:,:,ii) = imgaussfilt(backgroundStack(:,:,ii),gaussianSigma);
        end % for
    else
        gaussianStack = backgroundStack;
    end % if
    outName = strcat(outFolder,'03_GaussianFiltered_',tcrName);
    writeTiff(gaussianStack,outName);
    %% Determine +/- TCR signal
    structParameters.status = 'Status: Determine +/- TCR...';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
    means = zeros(dimensions(3),1);
    stds = zeros(dimensions(3),1);
    medians = zeros(dimensions(3),1);
    % If Auto is checked, use the contour mask
    if get(checkAuto,'Value') == 1
        for ii = 1:dimensions(3)
            temp = zeros(dimensions(1),dimensions(2));
            for jj = 1:dimensions(1)
                for kk = 1:dimensions(2)
                    if contourStack(jj,kk,ii) == 1
                        temp(jj,kk) = gaussianStack(jj,kk,ii);
                    else
                        temp(jj,kk) = NaN;
                    end % if
                end %
            end % for
            means(ii,1) = nanmean(nanmean(temp));
            medians(ii,1) = nanmedian(nanmedian(temp));
            stds(ii,1) = nanstd(nanstd(temp));
        end % for
    end % if
    % If Manual is checked, use the user input
    if get(checkManual,'Value') == 1
        for ii = 1:dimensions(3)
            temp = zeros((yMax-yMin+1),(xMax-xMin+1));
            for jj = yMin:yMax
                for kk = xMin:xMax
                    temp(jj-yMin+1,kk-xMin+1) = gaussianStack(jj,kk,ii);
                end % for
            end % for
            means(ii,1) = nanmean(nanmean(temp));
            medians(ii,1) = nanmedian(nanmedian(temp));
            stds(ii,1) = nanstd(nanstd(temp));
        end % for
    end % if
    % Determine cutoff
    if get(checkMean,'Value') == 1
        cutoffs = means+SDoverNoise*stds;
    elseif get(checkMedian,'Value') == 1
        cutoffs = medians+SDoverNoise*stds;
    end % if
    % Create threshold stacks
    threshIntensity = zeros(dimensions);
    threshBinary = zeros(dimensions);
    for ii = 1:dimensions(3)
        for jj = 1:dimensions(1)
            for kk = 1:dimensions(2)
                if gaussianStack(jj,kk,ii) > cutoffs(ii,1)
                    threshIntensity(jj,kk,ii) = gaussianStack(jj,kk,ii);
                    threshBinary(jj,kk,ii) = 1;
                end % if
            end % for
        end % for
    end % for
    outName = strcat(outFolder,'04_ThresholdIntensity_',tcrName);
    writeTiff(threshIntensity,outName);
    outName = strcat(outFolder,'05_BinaryIntensity_',tcrName);
    writeTiff(threshBinary,outName);
    % Remove small noise
    minusSpots = bwareaopen(threshBinary,3);
    minusSpots = uint8(minusSpots);
    outName = strcat(outFolder,'06_minusSpots_',tcrName);
    writeTiff(minusSpots,outName);
    %% Segmenting the contacts
    structParameters.status = 'Status: Separating Nanocontacts...';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
    % TCR positive nanocontacts
    TCRposNanocontacts = zeros(dimensions);
    NC = {};
    for ii = 1:dimensions(3)
        overlap = zeros(dimensions(1),dimensions(2));
        Ibw = watershedMaskStack(:,:,ii);
        CC = bwconncomp(Ibw,4);
        sample = threshBinary(:,:,ii);
        verts = reshape(sample,[dimensions(1)*dimensions(2),1]);
        for gg = 1:length(verts)
            if verts(gg,1) == 1
                for hh = 1:CC.NumObjects
                    for jj = 1:length(CC.PixelIdxList{1,hh})
                        if CC.PixelIdxList{1,hh}(jj,1) == gg
                            overlap(vertcat(CC.PixelIdxList{1,hh})) = true;
                            TCRposNanocontacts(:,:,ii) = overlap;
                        end % if
                    end % for
                end % for
            end % if
        end %
        NC.frames(ii) = CC;
        disp(ii)
    end % for
    
    outName = strcat(outFolder,'07_TCRpos_',tcrName);
    writeTiff(TCRposNanocontacts,outName);
    % TCR negative nanocontacts
    TCRnegNanocontacts = watershedMaskStack - TCRposNanocontacts;
    outName = strcat(outFolder,'08_TCRneg_',tcrName);
    writeTiff(TCRnegNanocontacts,outName);
    
    %% Box Method
    % TCR+
    sumStruct = zeros(dimensions);
    for jj = 1:(dimensions(1) - 5)
        for kk = 1:(dimensions(2) - 5)
            for ii = 1:dimensions(3)
                sumStruct(2+jj,2+kk,ii) = sum(sum(TCRposNanocontacts(jj:jj+4,kk:kk+4,ii)));
            end % for
        end % for
    end % for
    
    % TCR-
    sumNegStruct = zeros(dimensions);
    for jj = 1:(dimensions(1)-5)
        for kk = 1:(dimensions(2) - 5)
            for ii = 1:dimensions(3)
                sumNegStruct(2+jj,2+kk,ii) = sum(sum(TCRnegNanocontacts(jj:jj+4,kk:kk+4,ii)));
            end % for
        end % for
    end % for
    
    area = 13;
    % Record the length of all possible contacts
    % TCR+
    sumStruct = flip(sumStruct,3);
    contacts = zeros(dimensions);
    for jj = 1:dimensions(1)
        for kk = 1:dimensions(2)
            for ii = 1:dimensions(3)
                if sumStruct(jj,kk,ii) >= area
                    count = 1;
                    while (ii + count) <= dimensions(3) && sumStruct(jj,kk,ii+count) >= area
                        count = count + 1;
                    end % while
                    contacts(jj,kk,ii) = count;
                end % if
            end % for
        end % for
    end % for
    % TCR-
    sumNegStruct = flip(sumNegStruct,3);
    negContacts = zeros(dimensions);
    for jj = 1:dimensions(1)
        for kk = 1:dimensions(2)
            for ii = 1:dimensions(3)
                if sumNegStruct(jj,kk,ii) >= area
                    negCount = 1;
                    while (ii + negCount) <= dimensions(3) && sumNegStruct(jj,kk,ii+negCount) >= area
                        negCount = negCount + 1;
                    end % while
                    negContacts(jj,kk,ii) = negCount;
                end % if
            end % for
        end % for
    end % for
    % Removes pixels that are considered sequential contacts.
    % TCR+
    contactCount = zeros(dimensions);
    for jj = 1:dimensions(1)
        for kk = 1:dimensions(2)
            for ii = 1:dimensions(3)
                if contacts(jj,kk,ii) == 0
                    contactCount(jj,kk,ii) = 0;
                elseif (contacts(jj,kk,ii) > 0 && (ii + 1 <= dimensions(3) && contacts(jj,kk,ii+1) < contacts(jj,kk,ii)))
                    contactCount(jj,kk,ii) = contacts(jj,kk,ii);
                end % if
                if contacts(jj,kk,ii) > 0 && (ii-1 > 0 && contacts(jj,kk,ii) < contacts(jj,kk,ii-1))
                    contactCount(jj,kk,ii) = 0;
                end % if
            end % for
        end % for
    end % for
    % TCR-
    negContactCount = zeros(dimensions);
    for jj = 1:dimensions(1)
        for kk = 1:dimensions(2)
            for ii = 1:dimensions(3)
                if negContacts(jj,kk,ii) == 0
                    negContactCount(jj,kk,ii) = 0;
                elseif (negContacts(jj,kk,ii) > 0 && (ii + 1 <= dimensions(3) && negContacts(jj,kk,ii+1) < negContacts(jj,kk,ii)))
                    negContactCount(jj,kk,ii) = negContacts(jj,kk,ii);
                end % if
                if negContacts(jj,kk,ii) > 0 && (ii-1 > 0 && negContacts(jj,kk,ii) < negContacts(jj,kk,ii-1))
                    negContactCount(jj,kk,ii) = 0;
                end % if
            end % for
        end % for
    end % for
    % Reorient the stacks and save to .tif
    sumStruct = flip(sumStruct,3);
    outName = strcat(outFolder,'09_TCRposBoxMethod_',tcrName);
    writeTiff(sumStruct,outName); 
    sumNegStruct = flip(sumNegStruct,3);
    outName = strcat(outFolder,'10_TCRnegBoxMethod_',tcrName);
    writeTiff(sumNegStruct,outName);
    contacts = flip(contacts,3);
    outName = strcat(outFolder,'11_TCRposContactLength_',tcrName);
    writeTiff(contacts,outName);
    negContacts = flip(negContacts,3);
    outName = strcat(outFolder,'12_TCRnegContactLength_',tcrName);
    writeTiff(negContacts,outName);
    contactCount = flip(contactCount,3);
    outName = strcat(outFolder,'13_TCRposContactCount_',tcrName);
    writeTiff(contactCount,outName);
    negContactCount = flip(negContactCount,3);
    outName = strcat(outFolder,'14_TCRnegContactCount_',tcrName);
    writeTiff(negContactCount,outName);
    
    %% Get the TCR Intensity for the nanocontact pixels
    structParameters.status = 'Finding the Intensity of the Nanocontacts...';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
    % TCR+ and TCR-
    for ii = 1:dimensions(3)
        tempFrame = TCRStack(:,:,ii);
        tempLength = contacts(:,:,ii);
        tempNegLength = negContacts(:,:,ii);
        for hh = 1:NC.frames(ii).NumObjects
            for gg = 1:length(NC.frames(ii).PixelIdxList{1,hh})
                NC.frames(ii).TCRIntList{1,hh}(gg,1) = tempFrame(NC.frames(ii).PixelIdxList{1,hh}(gg,1));
                NC.frames(ii).ContactLengthList{1,hh}(gg,1) = tempLength(NC.frames(ii).PixelIdxList{1,hh}(gg,1));
                NC.frames(ii).negContactLengthList{1,hh}(gg,1) = tempNegLength(NC.frames(ii).PixelIdxList{1,hh}(gg,1));
            end % for
            NC.frames(ii).TCRmeans{1,hh} = mean(NC.frames(ii).TCRIntList{1,hh});
            NC.frames(ii).ContactLenFrames{1,hh} = max(NC.frames(ii).ContactLengthList{1,hh});
            NC.frames(ii).negContactLenFrames{1,hh} = max(NC.frames(ii).negContactLengthList{1,hh});
            NC.frames(ii).ContactLenSeconds{1,hh} = (max(NC.frames(ii).ContactLengthList{1,hh}))*timeRes;
            NC.frames(ii).negContactLenSeconds{1,hh} = (max(NC.frames(ii).negContactLengthList{1,hh}))*timeRes;
        end % for
    end % for        
    
    % Create tif stacks of the coded nanocontacts
    TCRmeans = zeros(dimensions(1),dimensions(2));
    contactPersistence = zeros(dimensions(1),dimensions(2));
    negContactPersistence = zeros(dimensions(1),dimensions(2));
    for ii = 1:dimensions(3)
        intFrame = zeros(1,dimensions(1)*dimensions(2));
        lenFrame = zeros(1,dimensions(1)*dimensions(2));
        lenNegFrame = zeros(1,dimensions(1)*dimensions(2));
        for hh = 1:NC.frames(ii).NumObjects
            for gg = 1:length(NC.frames(ii).PixelIdxList{1,hh})
                intFrame(NC.frames(ii).PixelIdxList{1,hh}(gg,1)) = NC.frames(ii).TCRmeans{1,hh};
                lenFrame(NC.frames(ii).PixelIdxList{1,hh}(gg,1)) = NC.frames(ii).ContactLenSeconds{1,hh};
                lenNegFrame(NC.frames(ii).PixelIdxList{1,hh}(gg,1)) = NC.frames(ii).negContactLenSeconds{1,hh};
            end % for
        end % for
        intFrame = reshape(intFrame,[dimensions(1),dimensions(2)]);
        lenFrame = reshape(lenFrame,[dimensions(1),dimensions(2)]);
        lenNegFrame = reshape(lenNegFrame,[dimensions(1),dimensions(2)]);
        TCRmeans(:,:,ii) = intFrame;
        contactPersistence(:,:,ii) = lenFrame;
        negContactPersistence(:,:,ii) = lenNegFrame;
    end % for
    
    outName = strcat(outFolder,'15_AllTCRmeans_',tcrName);
    writeTiff(TCRmeans,outName);
    outName = strcat(outFolder,'16_TCRposContactPersistence_',tcrName);
    writeTiff(contactPersistence,outName);
    outName = strcat(outFolder,'17_TCRnegContactPersistence_',tcrName);
    writeTiff(negContactPersistence,outName);
    
    TCRposList = [];
    TCRnegList = [];
    posCount = 0;
    negCount = 0;
    for ii = 1:dimensions(3)
        for hh = 1:NC.frames(ii).NumObjects
            if sum(NC.frames(ii).ContactLengthList{1,hh}) >= 1
                posCount = posCount + 1;
                TCRposList(posCount,1) = posCount;
                TCRposList(posCount,2) = NC.frames(ii).TCRmeans{1,hh};
                TCRposList(posCount,3) = NC.frames(ii).ContactLenSeconds{1,hh};
                TCRposList(posCount,4) = NC.frames(ii).ContactLenFrames{1,hh};
            elseif sum(NC.frames(ii).negContactLengthList{1,hh}) >= 1 
                negCount = negCount + 1;
                TCRnegList(negCount,1) = negCount;
                TCRnegList(negCount,2) = NC.frames(ii).TCRmeans{1,hh};
                TCRnegList(negCount,3) = NC.frames(ii).negContactLenSeconds{1,hh};
                TCRnegList(negCount,4) = NC.frames(ii).negContactLenFrames{1,hh};
            end % if
        end % for
    end % for
                
    %% Plot the TCR intensity vs dwell time
    
    hf2 = figure(2);
    scatter(TCRposList(:,3),TCRposList(:,2),[],'MarkerEdgeColor','r');
    xlabel('Dwell Time (s)');
    xlim([0 max(TCRposList(:,3))]);
    ylim([0 max(TCRposList(:,2))]);
    ylabel('TCR+ Intensity (mean)');
    
    hf3 = figure(3);
    scatter(TCRnegList(:,3),TCRnegList(:,2),[],'MarkerEdgeColor','r');
    xlabel('Dwell Time (s)');
    xlim([0 max(TCRposList(:,3))]);
    ylim([0 max(TCRposList(:,2))]);
    ylabel('TCR- Intensity (mean)');
    
    % Save some stats to the structure
    % TCR+
    NC.Summary.TCRpos.TotalContacts = length(TCRposList);
    NC.Summary.TCRpos.AveIntensity = mean(TCRposList(:,2));
    NC.Summary.TCRpos.AveDwellTimeS = mean(TCRposList(:,3));
    NC.Summary.TCRpos.AveDwellTimeFrames = mean(TCRposList(:,4));
    NC.Summary.TCRpos.StdDwellTimeS = std(TCRposList(:,3));
    NC.Summary.TCRpos.StdDwellTimeFrames = std(TCRposList(:,4));
    NC.Summary.TCRpos.MaxDwellTimeS = max(TCRposList(:,3));
    NC.Summary.TCRpos.MaxDwellTimeFrames = max(TCRposList(:,4));
     
    % TCR-
    NC.Summary.TCRneg.TotalContacts = length(TCRnegList);
    NC.Summary.TCRneg.AveIntensity = mean(TCRnegList(:,2));
    NC.Summary.TCRneg.AveDwellTimeS = mean(TCRnegList(:,3));
    NC.Summary.TCRneg.AveDwellTimeFrames = mean(TCRnegList(:,4));
    NC.Summary.TCRneg.StdDwellTimeS = std(TCRnegList(:,3));
    NC.Summary.TCRneg.StdDwellTimeFrames = std(TCRnegList(:,4));
    NC.Summary.TCRneg.MaxDwellTimeS = max(TCRnegList(:,3));
    NC.Summary.TCRneg.MaxDwellTimeFrames = max(TCRnegList(:,4));
    
    %% Savings
    % Save the figs
    saveas(hf2,strcat(outFolder,'TCRposDwellTime',tcrName(1:end-4)),'png'); 
    saveas(hf2,strcat(outFolder,'TCRposDwellTime',tcrName(1:end-4)),'fig');   
    saveas(hf3,strcat(outFolder,'TCRnegDwellTime',tcrName(1:end-4)),'png'); 
    saveas(hf3,strcat(outFolder,'TCRnegDwellTime',tcrName(1:end-4)),'fig');  
    % Save a .txt of contacts
    % TCR+
    TimeVsIntpos = TCRposList.';
    fid = fopen(strcat(outFolder,'TCRposNanocontacts.txt'),'wt');
    fprintf(fid,'ID, Mean Int, Dwell Time (s), Dwell Time Frames\n');
    fprintf(fid,'%g, %g, %g, %g\n',TimeVsIntpos);
    fclose(fid);
    % TCR-
    TimeVsIntneg = TCRnegList.';
    fid = fopen(strcat(outFolder,'TCRnegNanocontacts.txt'),'wt');
    fprintf(fid,'ID, Mean Int, Dwell Time (s), Dwell Time Frames\n');
    fprintf(fid,'%g, %g, %g, %g\n',TimeVsIntneg);
    fclose(fid);
    
    % Save a .txt of summary
    % TCR+
    fid = fopen(strcat(outFolder,'TCRposSummary.txt'),'wt');
    fprintf(fid,'Average Dwell Time (s), Std Dev Dwell Time (s), Max Dwell Time (s), Total Contacts\n');
    saveTxt = [NC.Summary.TCRpos.AveDwellTimeS, NC.Summary.TCRpos.StdDwellTimeS,NC.Summary.TCRpos.MaxDwellTimeS,NC.Summary.TCRpos.TotalContacts];
    fprintf(fid,'%g, %g, %g, %g\n',saveTxt);
    fprintf(fid,'Average Dwell Time (f), Std Dev Dwell Time (f), Max Dwell Time (f), Average Intensity\n');
    saveTxt2 = [NC.Summary.TCRpos.AveDwellTimeFrames, NC.Summary.TCRpos.StdDwellTimeFrames,NC.Summary.TCRpos.MaxDwellTimeFrames,NC.Summary.TCRpos.AveIntensity];
    fprintf(fid,'%g, %g, %g, %g\n',saveTxt2);
    fclose(fid);
    % TCR-
    fid = fopen(strcat(outFolder,'TCRnegSummary.txt'),'wt');
    fprintf(fid,'Average Dwell Time (s), Std Dev Dwell Time (s), Max Dwell Time (s), Total Contacts\n');
    saveTxt = [NC.Summary.TCRneg.AveDwellTimeS, NC.Summary.TCRneg.StdDwellTimeS,NC.Summary.TCRneg.MaxDwellTimeS,NC.Summary.TCRneg.TotalContacts];
    fprintf(fid,'%g, %g, %g, %g\n',saveTxt);
    fprintf(fid,'Average Dwell Time (f), Std Dev Dwell Time (f), Max Dwell Time (f), Average Intensity\n');
    saveTxt2 = [NC.Summary.TCRneg.AveDwellTimeFrames, NC.Summary.TCRneg.StdDwellTimeFrames,NC.Summary.TCRneg.MaxDwellTimeFrames,NC.Summary.TCRneg.AveIntensity];
    fprintf(fid,'%g, %g, %g, %g\n',saveTxt2);
    fclose(fid);
    % Save data in structures
    save(strcat(outFolder,'Nanocontacts.mat'),'NC'); 
    save(strcat(outFolder,'AnalysisParameters.mat'),'structParameters');
    
    structParameters.status = 'Status: Finished! Ready to run again!';
    setappdata(guiNC,'structParameters',structParameters);
    set(handles.textStatus,'String',structParameters.status);
    pause(0.01)
       
    toc
end % pushRun

