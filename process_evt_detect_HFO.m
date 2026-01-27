function varargout = process_evt_detect_HFO(varargin)
% PROCESS_EVT_DETECT_HFO: Detect High frequency oscillations in a continuous file, and create set of events called "HFO"

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
    sProcess.Comment     = 'Detect HFO in Spikes (within spike region)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 561;

    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.eventname.Comment = 'Event name:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'HFO_inSpikeRegion';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = cell(1, numel(sInputs));
    evtName = sProcess.options.eventname.Value;

    for f = 1:numel(sInputs)

        % ===== Load timefreq (Morlet) =====
        tfMat = in_bst(sInputs(f).FileName);
        morletTF = tfMat.TF;             % expected: [nChannels x nFreq x nTime]
        time = tfMat.Time(:)';           % seconds
        Fs = 1 / mean(diff(time));

        % ===== Load associated data epoch =====
        dataMat = in_bst(tfMat.DataFile);
        F = dataMat.F;                   % [nChannels x nSamples]
        time = dataMat.Time(:)';         % use data time

        nChannels = size(F,1);
        nSamples  = size(F,2);

        % ===== Matrices like in original code =====
        EoI = zeros(nChannels, nSamples);
        HFO = zeros(nChannels, nSamples);

        % ===== Block 2: Determine EoI =====
        for i = 1:nChannels
            channel = F(i,:);
            EoI(i,:) = determineEoI(channel, time);
        end

        % ===== Block 3: Determine HFO from EoI (needs TF map per channel) =====
        HFOrow = zeros(1, nSamples);
        for j = 1:nChannels
            time_frequency_map2 = squeeze(morletTF(j,:,:));   % [nFreq x nTime]
            EoIrow = EoI(j,:);
            channel = F(j,:);
            HFOrow = determineHFOfromEoI(channel, EoIrow, Fs, time, HFOrow, time_frequency_map2);
            HFO(j,:) = HFOrow;
            HFOrow = zeros(1, nSamples);
        end

        % ===== Block 4: GFP spike region (using F as "FNormal") =====
        [~, minimum_index, ~, min_after_spike] = computeGFPandSpikeRegion(F, time);

        % ===== Determine HFO in spike region and create events =====
        detTimes = [];
        detCh    = [];
        for i = 1:nChannels
            [~, yes_HFO, start_HFO, end_HFO] = determineHFOinSpikeRegion(HFO(i,:), time, minimum_index, min_after_spike);
            if yes_HFO == 1 && start_HFO ~= 0 && end_HFO ~= 0
                detTimes(end+1,:) = [start_HFO, end_HFO]; %#ok<AGROW>
                detCh(end+1,1)    = i;                  %#ok<AGROW>
            end
        end

        % ===== Add events to dataMat (minimal) =====
        if ~isfield(dataMat, 'Events') || isempty(dataMat.Events)
            dataMat.Events = repmat(struct('label',[],'color',[],'epochs',[],'times',[],'samples',[],'reactTimes',[],'select',[],'channels',[],'notes',[]), 0);
        end

        % find/create event
        iEvt = [];
        for k = 1:numel(dataMat.Events)
            if strcmpi(dataMat.Events(k).label, evtName)
                iEvt = k;
                break;
            end
        end
        if isempty(iEvt)
            iEvt = numel(dataMat.Events) + 1;
            dataMat.Events(iEvt).label  = evtName;
            dataMat.Events(iEvt).epochs = [];
            dataMat.Events(iEvt).times  = [];
            dataMat.Events(iEvt).samples = [];
            dataMat.Events(iEvt).channels = [];
            dataMat.Events(iEvt).select = 1;
        end

        if ~isempty(detTimes)
            t = detTimes';                       % [2 x nEvents]
            dataMat.Events(iEvt).epochs   = [dataMat.Events(iEvt).epochs, ones(1,size(detTimes,1))];
            dataMat.Events(iEvt).times    = [dataMat.Events(iEvt).times, t];
            dataMat.Events(iEvt).samples  = [dataMat.Events(iEvt).samples, round(t .* Fs)];
            dataMat.Events(iEvt).channels = [dataMat.Events(iEvt).channels, detCh(:)'];
        end

        % Save updated data epoch
        bst_save(file_fullpath(tfMat.DataFile), dataMat, 'v6');

        % Output: unchanged timefreq file (pipeline-friendly)
        OutputFiles{f} = sInputs(f).FileName;
    end
end

