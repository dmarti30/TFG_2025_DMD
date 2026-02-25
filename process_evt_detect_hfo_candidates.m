function varargout = process_evt_detect_hfo_candidates( varargin )
% PROCESS_EVT_DETECT_ECG: Detect heartbeats in a continuous file, and create set of events called "cardiac"

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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Detect HFO candidates';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 45;

    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Event name
    sProcess.options.eventname.Comment = 'Event name:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'HFO_candidate';

    % Channel name
    sProcess.options.channelname.Comment = 'Channel name:';
    sProcess.options.channelname.Type    = 'channelname';
    sProcess.options.channelname.Value   = '';

    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];

    % Robust threshold: thr_on = median(env) + k * 1.4826*MAD(env)
    sProcess.options.kmad.Comment = 'k (thr_on = median + k*MAD):';
    sProcess.options.kmad.Type    = 'value';
    sProcess.options.kmad.Value   = {5, '', 2};

    % Hysteresis: thr_off = median + (kOffRatio*k)*MAD
    sProcess.options.koff.Comment = 'Off ratio (thr_off = median + ratio*k*MAD):';
    sProcess.options.koff.Type    = 'value';
    sProcess.options.koff.Value   = {0.5, '', 2};

    % Minimum duration (ms)
    sProcess.options.mindur.Comment = 'Min duration (ms):';
    sProcess.options.mindur.Type    = 'value';
    sProcess.options.mindur.Value   = {20, 'ms', 0};

    % Merge gap (ms)
    sProcess.options.mergegap.Comment = 'Merge gap (ms):';
    sProcess.options.mergegap.Type    = 'value';
    sProcess.options.mergegap.Value   = {10, 'ms', 0};

    % Minimum cycles (approx) in the filtered waveform
    sProcess.options.mincycles.Comment = 'Min cycles:';
    sProcess.options.mincycles.Type    = 'value';
    sProcess.options.mincycles.Value   = {4, '', 0};

    % High-cut used to compute dynamic max zero-crossings (assume filter ~80-240)
    sProcess.options.fhigh.Comment = 'Assumed high-cut (Hz) for max-ZC rule:';
    sProcess.options.fhigh.Type    = 'value';
    sProcess.options.fhigh.Value   = {240, 'Hz', 0};

    % Max ZC tolerance factor (dynamic maxZC = 2*fhigh*dur*factor)
    sProcess.options.maxzcfactor.Comment = 'Max ZC factor (dynamic):';
    sProcess.options.maxzcfactor.Type    = 'value';
    sProcess.options.maxzcfactor.Value   = {1.3, '', 2};

    % Overwrite existing events with same name
    sProcess.options.overwrite.Comment = 'Overwrite existing events with same name';
    sProcess.options.overwrite.Type    = 'checkbox';
    sProcess.options.overwrite.Value   = 1;
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    evtName  = strtrim(sProcess.options.eventname.Value);
    chanName = strtrim(sProcess.options.channelname.Value);
    if isempty(evtName),  evtName  = 'HFO_candidate'; end
    if isempty(chanName)
        bst_report('Error', sProcess, [], 'Please select one channel.');
        return;
    end

    kMad       = sProcess.options.kmad.Value{1};
    kOffRatio  = sProcess.options.koff.Value{1};
    minDurMs   = sProcess.options.mindur.Value{1};
    mergeMs    = sProcess.options.mergegap.Value{1};
    minCycles  = sProcess.options.mincycles.Value{1};
    fHighAss   = sProcess.options.fhigh.Value{1};
    maxZCFact  = sProcess.options.maxzcfactor.Value{1};
    isOverwrite = sProcess.options.overwrite.Value;

    if isfield(sProcess.options,'timewindow') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    else
        TimeWindow = [];
    end

    % For each RAW file
    for iFile = 1:numel(sInputs)
        isRaw = strcmpi(sInputs(iFile).FileType, 'raw');
        if ~isRaw
            bst_report('Error', sProcess, sInputs(iFile), 'This process is RAW-only.');
            continue;
        end

        % === Load raw descriptor and time ===
        DataMat = in_bst_data(sInputs(iFile).FileName, 'F', 'Time');
        sFile   = DataMat.F;

        % Continuous only
        if ~isempty(sFile.epochs)
            bst_report('Error', sProcess, sInputs(iFile), 'This process only works on continuous RAW recordings (no epochs).');
            continue;
        end

        % Load channel file
        ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
        iChan = find(strcmpi({ChannelMat.Channel.Name}, chanName));
        if isempty(iChan)
            bst_report('Error', sProcess, sInputs(iFile), ['Channel "' chanName '" not found.']);
            continue;
        elseif numel(iChan) > 1
            bst_report('Error', sProcess, sInputs(iFile), ['Multiple channels named "' chanName '".']);
            continue;
        end

        % Time window bounds (samples in file)
        if ~isempty(TimeWindow)
            SamplesBounds = round(sFile.prop.times(1) .* sFile.prop.sfreq) + bst_closest(TimeWindow, DataMat.Time) - 1;
        else
            SamplesBounds = [];
        end

        % Read the selected channel from RAW
        [F, TimeVector] = in_fread(sFile, ChannelMat, 1, SamplesBounds, iChan);
        if isempty(F) || numel(TimeVector) < 2
            bst_report('Error', sProcess, sInputs(iFile), 'Invalid time window or empty read.');
            continue;
        end

        % Ensure row vector
        F = double(F(:))';
        TimeVector = double(TimeVector(:))';

        Fs = 1 / (TimeVector(2) - TimeVector(1));
        minDurSmp = max(1, round((minDurMs/1000) * Fs));
        mergeSmp  = max(0, round((mergeMs/1000)  * Fs));

        % === Hilbert envelope ===
        env = abs(hilbert(F));

        % === Robust threshold (median + k*MAD) ===
        m = median(env);
        madv = median(abs(env - m));
        sigmaRob = 1.4826 * madv;
        if sigmaRob <= eps
            % fallback if MAD collapses
            sigmaRob = std(env);
        end

        thr_on  = m + kMad * sigmaRob;
        thr_off = m + (kOffRatio * kMad) * sigmaRob;

        % === Find candidate segments using hysteresis ===
        [segStart, segEnd] = find_segments_hysteresis(env, thr_on, thr_off);

        % === Merge segments separated by small gaps ===
        [segStart, segEnd] = merge_segments_samples(segStart, segEnd, mergeSmp);

        % === Filter by min duration + zero-crossings min/max ===
        keep = false(1, numel(segStart));
        for k = 1:numel(segStart)
            s1 = segStart(k);
            s2 = segEnd(k);

            if (s2 - s1 + 1) < minDurSmp
                continue;
            end

            seg = F(s1:s2);
            zc = sum(diff(seg > 0) ~= 0);        % sign changes
            cycles = zc / 2;

            if cycles < minCycles
                continue;
            end

            dur = (s2 - s1 + 1) / Fs;
            maxZC = floor(2 * fHighAss * dur * maxZCFact);
            if maxZC < (2*minCycles)
                maxZC = 2*minCycles;
            end
            if zc > maxZC
                continue;
            end

            keep(k) = true;
        end

        segStart = segStart(keep);
        segEnd   = segEnd(keep);

        % === Convert to Brainstorm event times: [2 x nEvents] ===
        if isempty(segStart)
            bst_report('Info', sProcess, sInputs(iFile), ['No HFO candidates detected on "' chanName '".']);
            OutputFiles{end+1} = sInputs(iFile).FileName; %#ok<AGROW>
            continue;
        end

        evtTimes = [TimeVector(segStart); TimeVector(segEnd)];

        % === Create/overwrite event in RAW ===
        if ~isfield(sFile,'events') || isempty(sFile.events)
            sFile.events = repmat(db_template('event'), 0);
        end

        iEvt = find(strcmpi({sFile.events.label}, evtName), 1);

        if isempty(iEvt)
            iEvt = numel(sFile.events) + 1;
            sEvent = db_template('event');
            sEvent.label = evtName;
            sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
        else
            sEvent = sFile.events(iEvt);
            if isOverwrite
                sEvent.epochs = [];
                sEvent.times = [];
                sEvent.reactTimes = [];
                sEvent.channels = [];
                sEvent.notes = [];
            end
        end

        % Set event times (extended events: 2 x N)
        sEvent.times  = evtTimes;
        sEvent.epochs = ones(1, size(evtTimes,2));
        sEvent.channels = [];
        sEvent.notes = [];

        sFile.events(iEvt) = sEvent;

        % Save back to RAW file
        DataMat.F = sFile;
        DataMat = rmfield(DataMat, 'Time');
        bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);

        bst_report('Info', sProcess, sInputs(iFile), sprintf('Detected %d HFO candidates on %s.', size(evtTimes,2), chanName));

        OutputFiles{end+1} = sInputs(iFile).FileName; %#ok<AGROW>
    end
end


%% ===== LOCAL HELPERS =====
function [starts, ends] = find_segments_hysteresis(x, thr_on, thr_off)
    n = numel(x);
    starts = zeros(1,0);
    ends   = zeros(1,0);

    inSeg = false;
    s0 = 1;

    for i = 1:n
        if ~inSeg
            if x(i) > thr_on
                inSeg = true;
                s0 = i;
            end
        else
            if x(i) < thr_off
                inSeg = false;
                starts(end+1) = s0; %#ok<AGROW>
                ends(end+1)   = i-1; %#ok<AGROW>
            end
        end
    end

    if inSeg
        starts(end+1) = s0; %#ok<AGROW>
        ends(end+1)   = n;  %#ok<AGROW>
    end
end

function [starts, ends] = merge_segments_samples(starts, ends, mergeGapSmp)
    if isempty(starts)
        return;
    end
    sOut = starts(1);
    eOut = ends(1);

    for k = 2:numel(starts)
        if (starts(k) - eOut(end)) <= mergeGapSmp
            eOut(end) = max(eOut(end), ends(k));
        else
            sOut(end+1) = starts(k); %#ok<AGROW>
            eOut(end+1) = ends(k);   %#ok<AGROW>
        end
    end

    starts = sOut;
    ends   = eOut;
end