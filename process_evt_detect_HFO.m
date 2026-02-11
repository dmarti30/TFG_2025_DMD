function varargout = process_evt_detect_HFO(varargin)
% PROCESS_EVT_DETECT_HFO: Detect High frequency oscillations in spike epochs and create events.

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
    sProcess.Comment     = 'Detect HFO in spikes (spike region)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 561;

    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.eventname.Comment = 'Event name:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'HFO';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = cell(1, numel(sInputs));
    evtName = strtrim(sProcess.options.eventname.Value);
    if isempty(evtName)
        evtName = 'HFO';
    end

    % Keep original MATLAB detection functions unchanged: just make sure they are available.
    addCoreFunctionsToPath();

    for iFile = 1:numel(sInputs)
        % ===== Load Morlet TF file linked to the spike epoch =====
        tfMat = in_bst(sInputs(iFile).FileName);
        dataMat = in_bst(tfMat.DataFile);

        if ~isfield(tfMat, 'TF') || isempty(tfMat.TF)
            bst_report('Warning', sProcess, sInputs(iFile), 'Skipping file: TF matrix is empty.');
            OutputFiles{iFile} = sInputs(iFile).FileName;
            continue;
        end
        if ~isfield(dataMat, 'F') || isempty(dataMat.F)
            bst_report('Warning', sProcess, sInputs(iFile), 'Skipping file: data matrix is empty.');
            OutputFiles{iFile} = sInputs(iFile).FileName;
            continue;
        end

        % Expected TF layout: [nChannels x nFreq x nTime]
        morletTF = tfMat.TF;
        F = dataMat.F;
        timeData = dataMat.Time(:)';
        Fs = 1 / mean(diff(timeData));

        nChannels = size(F, 1);
        nSamples  = size(F, 2);

        % ===== Block 2: Determine EoI =====
        EoI = zeros(nChannels, nSamples);
        for iCh = 1:nChannels
            EoI(iCh, :) = determineEoI(F(iCh, :), timeData);
        end

        % ===== Block 3: Determine HFO from EoI =====
        HFO = zeros(nChannels, nSamples);
        for iCh = 1:nChannels
            tfMap = squeeze(morletTF(iCh, :, :));
            tfMap = abs(tfMap);
            tfMap = matchTfToDataLength(tfMap, nSamples);

            HFOrow = zeros(1, nSamples);
            HFO(iCh, :) = determineHFOfromEoI(F(iCh, :), EoI(iCh, :), Fs, timeData, HFOrow, tfMap);
        end

        % ===== Block 4: Compute GFP and spike region =====
        [~, minimum_index, ~, min_after_spike] = computeGFPandSpikeRegion(F, timeData);

        % ===== Keep only HFOs inside the spike region =====
        detTimes = zeros(0, 2);
        detCh = zeros(0, 1);

        for iCh = 1:nChannels
            [~, yes_HFO, start_HFO, end_HFO] = determineHFOinSpikeRegion(HFO(iCh, :), timeData, minimum_index, min_after_spike);
            if (yes_HFO == 1) && (start_HFO ~= 0) && (end_HFO ~= 0)
                detTimes(end+1, :) = [start_HFO, end_HFO]; %#ok<AGROW>
                detCh(end+1, 1) = iCh; %#ok<AGROW>
            end
        end

        % ===== Add/append Brainstorm event =====
        if ~isempty(detTimes)
            dataMat = appendExtendedEvent(dataMat, evtName, detTimes, detCh, Fs);
            bst_save(file_fullpath(tfMat.DataFile), dataMat, 'v6');
        end

        % Output is unchanged TF file (pipeline friendly)
        OutputFiles{iFile} = sInputs(iFile).FileName;
    end
end


%% ===== LOCAL HELPERS =====
function addCoreFunctionsToPath()
    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    coreDir = fullfile(thisDir, 'epilepsy_code');
    if exist(coreDir, 'dir') == 7
        addpath(coreDir);
    end
end

function tfMapOut = matchTfToDataLength(tfMapIn, nSamples)
    nTfTime = size(tfMapIn, 2);
    if nTfTime == nSamples
        tfMapOut = tfMapIn;
    elseif nTfTime > nSamples
        tfMapOut = tfMapIn(:, 1:nSamples);
    else
        tfMapOut = [tfMapIn, repmat(tfMapIn(:, end), 1, nSamples - nTfTime)];
    end
end

function dataMat = appendExtendedEvent(dataMat, evtName, detTimes, detCh, Fs)
    if ~isfield(dataMat, 'Events') || isempty(dataMat.Events)
        dataMat.Events = repmat(struct( ...
            'label', [], 'color', [], 'epochs', [], 'times', [], 'samples', [], ...
            'reactTimes', [], 'select', [], 'channels', [], 'notes', []), 0);
    end

    iEvt = [];
    for i = 1:numel(dataMat.Events)
        if strcmpi(dataMat.Events(i).label, evtName)
            iEvt = i;
            break;
        end
    end

    if isempty(iEvt)
        iEvt = numel(dataMat.Events) + 1;
        dataMat.Events(iEvt).label = evtName;
        dataMat.Events(iEvt).color = [1, 0, 0];
        dataMat.Events(iEvt).epochs = [];
        dataMat.Events(iEvt).times = [];
        dataMat.Events(iEvt).samples = [];
        dataMat.Events(iEvt).reactTimes = [];
        dataMat.Events(iEvt).select = 1;
        dataMat.Events(iEvt).channels = [];
        dataMat.Events(iEvt).notes = [];
    end

    t = detTimes';
    nDet = size(detTimes, 1);

    dataMat.Events(iEvt).epochs = [dataMat.Events(iEvt).epochs, ones(1, nDet)];
    dataMat.Events(iEvt).times = [dataMat.Events(iEvt).times, t];
    dataMat.Events(iEvt).samples = [dataMat.Events(iEvt).samples, round(t .* Fs)];

    % Keep channels as numeric indices (one per detected interval), consistent with detection output.
    dataMat.Events(iEvt).channels = [dataMat.Events(iEvt).channels, detCh(:)'];
end
