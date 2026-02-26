function varargout =process_evt_detect_hfo_candidates_by_channel( varargin )
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
    sProcess.Comment     = 'Detect HFO candidates by channel';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 47;

    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Event name prefix
    sProcess.options.prefix.Comment = 'Event label prefix:';
    sProcess.options.prefix.Type    = 'text';
    sProcess.options.prefix.Value   = 'HFOcand';

    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];

    % Threshold: thr_on = median(env) + k * 1.4826*MAD(env)
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

    % High-cut used to compute dynamic max zero-crossings
    sProcess.options.fhigh.Comment = 'Assumed high-cut (Hz) for max-ZC rule:';
    sProcess.options.fhigh.Type    = 'value';
    sProcess.options.fhigh.Value   = {240, 'Hz', 0};

    % Max ZC tolerance factor (dynamic maxZC = 2*fhigh*dur*factor)
    sProcess.options.maxzcfactor.Comment = 'Max ZC factor (dynamic):';
    sProcess.options.maxzcfactor.Type    = 'value';
    sProcess.options.maxzcfactor.Value   = {1.3, '', 2};

    % Overwrite existing events with same per-channel labels
    sProcess.options.overwrite.Comment = 'Overwrite existing per-channel events';
    sProcess.options.overwrite.Type    = 'checkbox';
    sProcess.options.overwrite.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = cell(1, numel(sInputs));

    prefix     = strtrim(sProcess.options.prefix.Value);
    if isempty(prefix), prefix = 'HFOcand'; end

    kMad       = sProcess.options.kmad.Value{1};
    kOffRatio  = sProcess.options.koff.Value{1};
    minDurMs   = sProcess.options.mindur.Value{1};
    mergeMs    = sProcess.options.mergegap.Value{1};
    minCycles  = sProcess.options.mincycles.Value{1};
    fHighAss   = sProcess.options.fhigh.Value{1};
    maxZCFact  = sProcess.options.maxzcfactor.Value{1};
    isOverwrite = sProcess.options.overwrite.Value;

    % Time window
    TimeWindow = [];
    if isfield(sProcess.options,'timewindow') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    end

    for iFile = 1:numel(sInputs)
        try
            % Load RAW descriptor (sFile)
            DataMat = in_bst_data(sInputs(iFile).FileName, 'F');
            sFile   = DataMat.F;

            % Continuous RAW only
            if ~isempty(sFile.epochs)
                bst_report('Error', sProcess, sInputs(iFile), 'RAW must be continuous (no epochs).');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            % Load channel file
            ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
            chanNames  = {ChannelMat.Channel.Name};
            chanTypes  = {ChannelMat.Channel.Type};

            % Prefer EEG channels; fallback to all channels if none tagged EEG
            iChanList = find(strcmpi(chanTypes, 'EEG'));
            if isempty(iChanList)
                iChanList = 1:numel(ChannelMat.Channel);
            end

            % Read all selected channels from RAW
            [F, Time] = in_fread(sFile, ChannelMat, 1, [], iChanList);
            if isempty(F) || numel(Time) < 2
                bst_report('Error', sProcess, sInputs(iFile), 'Could not read RAW data.');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            F    = double(F);        % [nChan x nTime]
            Time = double(Time(:)'); % [1 x nTime]
            Fs   = 1 / (Time(2) - Time(1));

            % Apply time window if provided
            if ~isempty(TimeWindow)
                i1 = bst_closest(TimeWindow(1), Time);
                i2 = bst_closest(TimeWindow(2), Time);
                if i2 < i1, tmp = i1; i1 = i2; i2 = tmp; end
                F = F(:, i1:i2);
                Time = Time(i1:i2);
                if numel(Time) < 2
                    bst_report('Error', sProcess, sInputs(iFile), 'Invalid time window (empty selection).');
                    OutputFiles{iFile} = sInputs(iFile).FileName;
                    continue;
                end
            end

            minDurSmp = max(1, round((minDurMs/1000) * Fs));
            mergeSmp  = max(0, round((mergeMs/1000)  * Fs));

            % Ensure events exist
            if ~isfield(sFile,'events') || isempty(sFile.events)
                sFile.events = repmat(db_template('event'), 0);
            end

            % Process each channel and create an event label per channel
            for ii = 1:size(F,1)
                x = F(ii,:);

                % Hilbert envelope
                env = abs(hilbert(x));

                % Robust threshold (median + k*MAD)
                m = median(env);
                madv = median(abs(env - m));
                sigmaRob = 1.4826 * madv;
                if sigmaRob <= eps
                    sigmaRob = std(env);
                end
                thr_on  = m + kMad * sigmaRob;
                thr_off = m + (kOffRatio * kMad) * sigmaRob;

                % Candidate segments using hysteresis
                [segStart, segEnd] = find_segments_hysteresis(env, thr_on, thr_off);

                % Merge close segments
                [segStart, segEnd] = merge_segments_samples(segStart, segEnd, mergeSmp);

                % Filter by min duration + ZC min/max
                keepStart = zeros(1,0);
                keepEnd   = zeros(1,0);

                for k = 1:numel(segStart)
                    s1 = segStart(k);
                    s2 = segEnd(k);

                    if (s2 - s1 + 1) < minDurSmp
                        continue;
                    end

                    seg = x(s1:s2);

                    % zero-crossings (sign changes)
                    zc = sum(diff(seg > 0) ~= 0);
                    cycles = zc / 2;

                    if cycles < minCycles
                        continue;
                    end

                    dur = (s2 - s1 + 1) / Fs;
                    maxZC = floor(2 * fHighAss * dur * maxZCFact);
                    if maxZC < 2*minCycles
                        maxZC = 2*minCycles;
                    end
                    if zc > maxZC
                        continue;
                    end

                    keepStart(end+1) = s1; %#ok<AGROW>
                    keepEnd(end+1)   = s2; %#ok<AGROW>
                end

                % Build per-channel event label
                chLabel = sanitize_label(chanNames{iChanList(ii)});
                evtLabel = [prefix '_' chLabel];

                % Find/create this per-channel event
                iEvt = find(strcmpi({sFile.events.label}, evtLabel), 1);

                % Overwrite existing per-channel event (optional)
                if ~isempty(iEvt) && isOverwrite
                    sFile.events(iEvt).times = [];
                    sFile.events(iEvt).epochs = [];
                    sFile.events(iEvt).reactTimes = [];
                    sFile.events(iEvt).channels = [];
                    sFile.events(iEvt).notes = [];
                end

                % If no detections for this channel, skip creating new event
                if isempty(keepStart)
                    continue;
                end

                evtTimes = [Time(keepStart); Time(keepEnd)]; % [2 x N]

                if isempty(iEvt)
                    iEvt = numel(sFile.events) + 1;
                    sEvent = db_template('event');
                    sEvent.label = evtLabel;
                    sEvent.color = panel_record('GetNewEventColor', iEvt, sFile.events);
                else
                    sEvent = sFile.events(iEvt);
                end

                % Assign events (extended)
                sEvent.times  = [sEvent.times, evtTimes];
                sEvent.epochs = [sEvent.epochs, ones(1, size(evtTimes,2))];
                sEvent.channels = [];
                sEvent.notes = [];

                sFile.events(iEvt) = sEvent;
            end

            % Save back to RAW file
            DataMat.F = sFile;
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);

            bst_report('Info', sProcess, sInputs(iFile), 'HFO candidates by channel: done.');
            OutputFiles{iFile} = sInputs(iFile).FileName;

        catch ME
            loc = '';
            if ~isempty(ME.stack)
                loc = sprintf(' (%s:%d)', ME.stack(1).name, ME.stack(1).line);
            end
            bst_report('Error', sProcess, sInputs(iFile), sprintf('Error: %s%s', ME.message, loc));
            OutputFiles{iFile} = sInputs(iFile).FileName;
        end
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

function s = sanitize_label(name)
    % Replace non-alphanumeric with underscore, avoid empty labels
    s = regexprep(name, '[^A-Za-z0-9]', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_+', '');
    if isempty(s)
        s = 'Chan';
    end
end
