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

    OutputFiles = cell(1, numel(sInputs));
    evtName = sProcess.options.eventname.Value;

    for f = 1:numel(sInputs)

        %% === LOAD TIMEFREQ (MORLET) ===
        tfMat = in_bst(sInputs(f).FileName);
        TF   = tfMat.TF;              % [channels x freq x time]
        time = tfMat.Time(:)';        % time vector

        %% === LOAD ASSOCIATED DATA (SPIKE EPOCH) ===
        dataMat = in_bst(tfMat.DataFile);
        F = dataMat.F;                % [channels x samples]
        time = dataMat.Time(:)';      % overwrite with data time
        Fs = 1 / mean(diff(time));

        nChannels = size(F,1);
        nSamples  = size(F,2);

        %% === INIT MATRICES (LIKE TFG SCRIPT) ===
        EoI = zeros(nChannels, nSamples);
        HFO = zeros(nChannels, nSamples);

        %% === BLOCK 2: DETERMINE EoI ===
        for ch = 1:nChannels
            EoI(ch,:) = determineEoI(F(ch,:), time);
        end

        %% === BLOCK 3: DETERMINE HFO FROM EoI ===
        for ch = 1:nChannels
            time_frequency_map = squeeze(TF(ch,:,:));   % [freq x time]
            HFOrow = zeros(1, nSamples);
            HFOrow = determineHFOfromEoI(F(ch,:), EoI(ch,:), Fs, time, HFOrow, time_frequency_map);
            HFO(ch,:) = HFOrow;
        end

        %% === BLOCK 4: COMPUTE GFP & SPIKE REGION ===
        [~, minimum_index, ~, min_after_spike] = computeGFPandSpikeRegion(F, time);

        %% === DETERMINE HFO WITHIN SPIKE REGION & CREATE EVENTS ===
        detTimes = [];
        detCh    = [];

        for ch = 1:nChannels
            [~, yes_HFO, start_HFO, end_HFO] = ...
                determineHFOinSpikeRegion(HFO(ch,:), time, minimum_index, min_after_spike);

            if yes_HFO == 1 && start_HFO ~= 0 && end_HFO ~= 0
                detTimes(end+1,:) = [start_HFO, end_HFO]; %#ok<AGROW>
                detCh(end+1)      = ch;                  %#ok<AGROW>
            end
        end

        %% === ADD EVENTS TO DATA ===
        if ~isempty(detTimes)

            if ~isfield(dataMat, 'Events') || isempty(dataMat.Events)
                dataMat.Events = struct([]);
            end

            evtIdx = [];
            for k = 1:numel(dataMat.Events)
                if strcmp(dataMat.Events(k).label, evtName)
                    evtIdx = k;
                end
            end

            if isempty(evtIdx)
                evtIdx = numel(dataMat.Events) + 1;
                dataMat.Events(evtIdx).label    = evtName;
                dataMat.Events(evtIdx).epochs   = [];
                dataMat.Events(evtIdx).times    = [];
                dataMat.Events(evtIdx).samples  = [];
                dataMat.Events(evtIdx).channels = [];
                dataMat.Events(evtIdx).select   = 1;
            end

            dataMat.Events(evtIdx).epochs   = [dataMat.Events(evtIdx).epochs, ones(1,size(detTimes,1))];
            dataMat.Events(evtIdx).times    = [dataMat.Events(evtIdx).times, detTimes'];
            dataMat.Events(evtIdx).samples  = [dataMat.Events(evtIdx).samples, round(detTimes' * Fs)];
            dataMat.Events(evtIdx).channels = [dataMat.Events(evtIdx).channels, detCh];

            bst_save(file_fullpath(tfMat.DataFile), dataMat, 'v6');
        end

        OutputFiles{f} = sInputs(f).FileName;
    end
end
