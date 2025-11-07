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
% Authors: Francois Tadel, Elizabeth Bock, 2011-2016

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
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>   
% Create a new event that occurs at 1 second for all channels
eventTime = 1; % Time in seconds
eventName = sProcess.options.eventname.Value; % Get the event name from options

% Create the event structure
newEvent = struct();
newEvent.label = eventName; % Set the event label
newEvent.time = eventTime; % Set the event time

% Initialize the output files
OutputFiles = [];

% Add the new event to the output
i = 1; % Initialize index for OutputFiles
OutputFiles(i).event = newEvent; % Assign the new event to each channel

% Ensure that the output files have a filename field
if ~isfield(OutputFiles, 'FileName')
    OutputFiles(i).FileName = ''; % Assign an empty string if FileName is not present
end
% Save the new event in the events structure
if ~isfield(OutputFiles, 'events')
    OutputFiles(i).events = []; % Initialize events field if it doesn't exist
end

% Append the new event to the events field
OutputFiles(i).events = [OutputFiles(i).events; newEvent];
end

