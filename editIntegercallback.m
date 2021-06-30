function editIntegercallback(hObject, ~, hObjectContainer, gui)
%editProtrusionSizecallback Updates edit boxes and ensures they are integers
% R2015b
% 
% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% Archived March 2017

    %% Get the parameters structure
    structParameters = getappdata(gui,'structParameters');
    
    %% Get the converted string for the calling editbox
    newValue = str2double(get(hObject,'String'));
    
    %% Test for a valid value
    % Input must be a number greater than zero
    if isnan(newValue) || newValue <= 0
        set(hObject, 'String', hObjectContainer.OldString)
    else 
        % Test to see if the number is an integer
        if ~mod(newValue,1) == 1
            % Update the control property
            hObjectContainer.OldString = newValue;
            
            % Update the appdata
            if strcmp(hObject.Tag,'editTophat') == 1
                structParameters.backgroundDisc = newValue;
                setappdata(gui,'structParameters',structParameters)
            elseif strcmp(hObject.Tag,'editGaussianSigma') == 1
                structParameters.gaussianSigma = newValue;
                setappdata(gui,'structParameters',structParameters)
            elseif strcmp(hObject.Tag,'editSD') == 1
                structParameters.SDoverNoise = newValue;
                setappdata(gui,'structParameters',structParameters)
            end
                
        else
            % If not a whole number, round
            newValue = round(newValue);
            
            % Update the control property
            hObjectContainer.OldString = newValue;
            
            % Update the appdata
            if strcmp(hObject.Tag,'editTophat') == 1
                structParameters.backgroundDisc = newValue;
                setappdata(gui,'structParameters',structParameters)
            elseif strcmp(hObject.Tag,'editGaussianSigma') == 1
                structParameters.gaussianSigma = newValue;
                setappdata(gui,'structParameters',structParameters)
            elseif strcmp(hObject.Tag,'editSD') == 1
                structParameters.SDoverNoise = newValue;
                setappdata(gui,'structParameters',structParameters)
            end
            % Reset the string visible in the GUI
            set(hObject, 'String', hObjectContainer.OldString)
            setappdata(gui,'structParameters',structParameters)
        end % if
    end % if
            

    
end % editIntegercallback

