function varargout = process_evt_detect_EoI( varargin )
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
    sProcess.Comment     = 'Detect EoI';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'HFO detection';
    sProcess.Index       = 44;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsDetect#Detection:_Heartbeats';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Channel name
    sProcess.options.channelname.Comment = 'Channel name: ';
    sProcess.options.channelname.Type    = 'channelname';
    sProcess.options.channelname.Value   = '';
    % Channel name comment
    sProcess.options.channelhelp.Comment = '<I><FONT color="#777777">You can use the montage syntax here: "ch1, -ch2"</FONT></I>';
    sProcess.options.channelhelp.Type    = 'label';
    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Event name
    sProcess.options.eventname.Comment = 'Event name: ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'EoI';
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
end

function[EoI] = determineEoI(channel,time, ~)
    % Compute the envelope
    y = hilbert(channel); % Apply Hilbert transform to the channel data
    env = abs(y); % Compute the envelope

    % Standard Deviation and Mean of the envelope
    env_mean = mean(env);
    env_SD = std(env);

    % Set the threshold
    thr = env_mean + 1 * env_SD;
    duration_thr = 0.5 * thr;

    % Mark The EoI
    EoI = zeros(size(env));
    for i = 1:length(env)
        if env(i) >= thr
            EoI(i) = channel(i);
        else
            EoI(i) = 0;
        end
    end

    % I don't want to study the EoI above the spike
    for i = 1:length(EoI)
        if time(i) > 0.1
            EoI(i) = 0;
        end
    end

    % Look for the duration of the EoI
    for i = 1:length(EoI)
        if EoI(i) ~= 0
            value = 0;
            rest = 1;
            while value == 0 && (i - rest) >= 1 
                if env(i - rest) >= duration_thr || env(i - rest) <= -duration_thr 
                    EoI(i - rest) = channel(i - rest);
                    rest = rest + 1;
                else
                    value = 1;
                end
            end
        end
    end

    for i = 1:length(EoI)
        if EoI(i) ~= 0
            value = 0;
            sum = 1;
            while value == 0 && (i + sum) <= length(EoI)
                if env(i + sum) >= duration_thr || env(i + sum) <= -duration_thr 
                    EoI(i + sum) = channel(i + sum);
                    sum = sum + 1;
                else
                    value = 1;
                end
            end
        end
    end

    % Merge together EoI in an interval less than 30ms (10 samples in between)
    interval_threshold = (abs(time(1)) - abs(time(10)));
    i = 1;
    while i <= length(EoI)
        if EoI(i) ~= 0 && i < length(EoI) 
            a = 0;
            z = i + 1;
            cont = i;
            while a == 0 && z <= length(EoI)
                if EoI(z) ~= 0
                    z = z + 1;
                    a = 0;
                else
                    a = 1;
                    i = z - 1;
                    cont = z - 1;
                end
            end
            b = cont + 1;
            c = 0;
            while b <= length(EoI) && c == 0
                if EoI(b) ~= 0 
                    if time(b) < 0 && time(cont) < 0
                        interval = abs(time(cont)) - abs(time(b));
                    elseif time(cont) < 0 && time(b) == 0
                        interval = abs(time(cont));
                    elseif time(cont) < 0 && time(b) > 0
                        interval = abs(time(cont)) + time(b);
                    elseif time(cont) > 0 && time(b) > 0
                        interval = time(b) - time(cont);
                    else
                        interval = 10000;
                    end
                    c = 1;
                    if interval <= interval_threshold
                        for j = (cont + 1):b
                            if EoI(j) == 0
                                EoI(j) = channel(j);
                            end
                        end
                    end
                end
                b = b + 1;
            end
        end
        i = i + 1;
    end
    
    % Create the event structure
    sEvent = db_template('event');
    sEvent.label = sProcess.options.eventname.Value; % Get the event name from options
    sEvent.color = [1 0 0];
    sEvent.epochs = []; % Initialize epochs field
    sEvent.times = time(EoI ~= 0); % Set the event times based on EoI
    sEvent.reactTimes = []; % Initialize reactTimes field
    sEvent.select = []; % Initialize select field
    sEvent.channels = []; % Specify channels if needed, empty for all channels
    sEvent.notes = []; % Initialize notes field

    % Append the new event to the events field
    sFile.events(end+1) = sEvent;

    % Save the EoI results to output
    OutputFiles{iFile} = bst_save(sFile, EoI, 'EoI', sProcess.options.eventname.Value);
end
end