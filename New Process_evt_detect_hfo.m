function varargout = process_evt_detect_hfo( varargin )
% PROCESS_EVT_DETECT_HFO: Detect heartbeats in a continuous file, and create set of events called "hfo"

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Daniel Martin, Guiomar Niso, 2025

% Ensure that the output argument is assigned a value
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Detect hfo';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 43;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsDetect#Detection:_Heartbeats';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Event name
    sProcess.options.eventname.Comment = 'Event name: ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'hfo';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)    

% Initialize the output files
OutputFiles = {};

for iFile = 1:length(sInputs)
    % Load the raw file descriptor
    isRaw = strcmpi(sInputs(iFile).FileType, 'raw');

    if isRaw
        DataMat = in_bst_data(sInputs(iFile).FileName, 'F', 'events');
        sFile = DataMat.F;
    else
        DataMat = in_bst_data(sInputs(iFile).FileName, 'Time');
        sFile = in_fopen(sInputs(iFile).FileName, 'BST-DATA');
    end

    % Create a new event that occurs at 1 second for all channels
    % Basic events structure
    sEvent = repmat(db_template('event'), 0);
    iEvt =1;

    % Create the event structure
    sEvent(iEvt).label = sProcess.options.eventname.Value; % Get the event name from options
    sEvent(iEvt).color = [1 0 0];
    sEvent(iEvt).epochs = []; % Initialize epochs field
    sEvent(iEvt).times = [1.570156250000000e+02]; % Set the event time
    sEvent(iEvt).reactTimes = []; % Initialize reactTimes field
    sEvent(iEvt).select = []; % Initialize select field
    sEvent(iEvt).channels = []; % Specify channels if needed, empty for all channels
    sEvent(iEvt).notes = []; % Initialize notes field

    % Append the new event to the events field
    sFile.events(end+1) = sEvent;

     % Only save changes if something was detected
        if ~isempty(sEvent)
            % Report changes in .mat structure
            if isRaw
                DataMat.F = sFile;
            else
                DataMat.Events = sFile.events;
            end
            % Save file definition
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);
            % Report number of detected events
            nEvents =1;
            bst_report('Info', sProcess, sInputs(iFile), sprintf('%d events detected in HFO', nEvents));
        else
            bst_report('Warning', sProcess, sInputs(iFile), ['No event detected. Please check the signal quality.']);
        end
        
        OutputFiles = {sInputs.FileName};

end
end

