function pushMaskFilecallback( ~,~,guiNC,handles )
%pushInputFilecallback Prompts the user to find the .tif file to be
%analyzed. R2015b
%
% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% Archived March 2017

    %% Prompt the user to select a file
    % check for a currently selected file
    structParameters = getappdata(guiNC,'structParameters');
    inputFile = structParameters.maskName;
    
    if ~isempty(inputFile)
        [fileSelection,pathSelection] = uigetfile({'*.tif'},'Select .tif',structParameters.pathDir);
    else
        [fileSelection,pathSelection] = uigetfile({'*.tif'},'Select .tif');
    end % if

    %% Check for a canceled selection.
    if fileSelection(1) == 0 || strcmp(fileSelection,inputFile)
        return
    end % if
    
    %% Store the file as appdata
    structParameters.maskName = fileSelection;
    setappdata(guiNC,'structParameters',structParameters);
    handles.displayMaskFile.String = structParameters.maskName;


end

