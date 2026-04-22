function varargout = process_evt_detect_hfos_original_raw_no_bads(varargin)
% PROCESS_EVT_DETECT_HFOS_FROM_RAW_TO_RAW:
% Starting from the ORIGINAL raw file:
%   1) Apply notch filter in memory: 50,100,150,200 Hz
%   2) Apply band-pass filter in memory: 60-250 Hz
%   3) Run EXACTLY the same HFO candidate detector on the filtered signal
%   4) Save the detected events on the ORIGINAL raw file
%   5) Optionally also save the intermediate filtered files in Brainstorm
%
% IMPORTANT:
%   - The HFO detection logic below is intentionally preserved from the
%     original process_evt_detect_hfos_without_bad.
%   - The only addition is the preprocessing stage before detection.

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Detect HFO from original RAW (notch+band, save on original)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 54;

    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Event name
    sProcess.options.eventname.Comment = 'Event name:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'HFO_candidate';

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

    % Minimum cycles in the filtered waveform
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

    % Overwrite existing events with same name
    sProcess.options.overwrite.Comment = 'Overwrite existing events with same name';
    sProcess.options.overwrite.Type    = 'checkbox';
    sProcess.options.overwrite.Value   = 1;

    % ===== NEW: preprocessing options =====
    sProcess.options.notchfreqs.Comment = 'Notch frequencies (Hz):';
    sProcess.options.notchfreqs.Type    = 'value';
    sProcess.options.notchfreqs.Value   = {[50 100 150 200], 'Hz', 0};

    sProcess.options.highpass.Comment = 'Band-pass high cutoff (Hz):';
    sProcess.options.highpass.Type    = 'value';
    sProcess.options.highpass.Value   = {60, 'Hz', 0};

    sProcess.options.lowpass.Comment = 'Band-pass low cutoff (Hz):';
    sProcess.options.lowpass.Type    = 'value';
    sProcess.options.lowpass.Value   = {250, 'Hz', 0};

    sProcess.options.keepfiltered.Comment = 'Keep filtered files in Brainstorm database';
    sProcess.options.keepfiltered.Type    = 'checkbox';
    sProcess.options.keepfiltered.Value   = 0;
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
        evtName = 'HFO_candidate';
    end

    kMad        = sProcess.options.kmad.Value{1};
    kOffRatio   = sProcess.options.koff.Value{1};
    minDurMs    = sProcess.options.mindur.Value{1};
    mergeMs     = sProcess.options.mergegap.Value{1};
    minCycles   = sProcess.options.mincycles.Value{1};
    fHighAss    = sProcess.options.fhigh.Value{1};
    maxZCFact   = sProcess.options.maxzcfactor.Value{1};
    isOverwrite = sProcess.options.overwrite.Value;

    NotchFreqs  = sProcess.options.notchfreqs.Value{1};
    HighPassHz  = sProcess.options.highpass.Value{1};
    LowPassHz   = sProcess.options.lowpass.Value{1};
    KeepFiltered = sProcess.options.keepfiltered.Value;

    % Time window
    TimeWindow = [];
    if isfield(sProcess.options,'timewindow') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    end

    for iFile = 1:numel(sInputs)
        try
            % ===== LOAD ORIGINAL RAW =====
            DataMat = in_bst_data(sInputs(iFile).FileName, 'F');
            sFile   = DataMat.F;

            if ~isempty(sFile.epochs)
                bst_report('Error', sProcess, sInputs(iFile), 'RAW must be continuous (no epochs).');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
            chanNames  = {ChannelMat.Channel.Name};
            chanTypes  = {ChannelMat.Channel.Type};

            % Same channel selection policy as original detector
            iChanList = find(strcmpi(chanTypes, 'EEG'));
            if isempty(iChanList)
                iChanList = 1:numel(ChannelMat.Channel);
            end

            % Read ORIGINAL raw
            [Fraw, Time] = in_fread(sFile, ChannelMat, 1, [], iChanList);
            if isempty(Fraw) || numel(Time) < 2
                bst_report('Error', sProcess, sInputs(iFile), 'Could not read RAW data.');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            Fraw = double(Fraw);         % [nChan x nTime]
            Time = double(Time(:)');     % [1 x nTime]
            Fs   = 1 / (Time(2) - Time(1));

            % Apply time window if provided
            if ~isempty(TimeWindow)
                i1 = bst_closest(TimeWindow(1), Time);
                i2 = bst_closest(TimeWindow(2), Time);
                if i2 < i1
                    tmp = i1; i1 = i2; i2 = tmp;
                end
                Fraw = Fraw(:, i1:i2);
                Time = Time(i1:i2);
                if numel(Time) < 2
                    bst_report('Error', sProcess, sInputs(iFile), 'Invalid time window (empty selection).');
                    OutputFiles{iFile} = sInputs(iFile).FileName;
                    continue;
                end
                Fs = 1 / (Time(2) - Time(1));
            end

            % ===== PREPROCESSING IN MEMORY =====
            % 1) Notch
            Ffilt = Fraw;
            if ~isempty(NotchFreqs)
                Ffilt = process_notch('Compute', Ffilt, Fs, NotchFreqs);
            end

            % 2) Band-pass
            % Supported external call documented by Brainstorm:
            % process_bandpass('Compute', x, Fs, HighPass, LowPass, 'bst-hfilter', isMirror, isRelax)
            Ffilt = process_bandpass('Compute', Ffilt, Fs, HighPassHz, LowPassHz, 'bst-hfilter', 0, 0);

            % ===== OPTIONAL: SAVE FILTERED FILES IN DATABASE =====
            % This does NOT affect the detector or where events are written.
            if KeepFiltered
                try
                    sRawIn = {sInputs(iFile).FileName};

                    sFilesNotch = bst_process('CallProcess', 'process_notch', sRawIn, [], ...
                        'freqlist',    NotchFreqs, ...
                        'sensortypes', 'EEG', ...
                        'read_all',    1, ...
                        'overwrite',   0);

                    bst_process('CallProcess', 'process_bandpass', sFilesNotch, [], ...
                        'sensortypes',  'EEG', ...
                        'highpass',     HighPassHz, ...
                        'lowpass',      LowPassHz, ...
                        'attenuation',  'strict', ...
                        'ver',          '2019', ...
                        'mirror',       0, ...
                        'read_all',     1, ...
                        'overwrite',    0);
                catch MEkeep
                    bst_report('Warning', sProcess, sInputs(iFile), ...
                        sprintf('Filtered files were not saved, but detection continues normally: %s', MEkeep.message));
                end
            end

            % ===== FROM HERE, DETECTION LOGIC IS KEPT AS IN THE ORIGINAL SCRIPT =====
            minDurSmp = max(1, round((minDurMs/1000) * Fs));
            mergeSmp  = max(0, round((mergeMs/1000)  * Fs));

            % BAD intervals are read from the ORIGINAL raw events
            badTimes = get_bad_intervals(sFile.events);

            % Collect detections across channels
            allStart = [];
            allEnd   = [];
            allNote  = {};

            % Detect per channel
            for ii = 1:size(Ffilt,1)
                x = Ffilt(ii,:);

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

                % Filter by min duration + ZC min/max + BAD overlap
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

                    % reject if overlapping any BAD segment
                    evtStart = Time(s1);
                    evtEnd   = Time(s2);
                    if overlaps_bad(evtStart, evtEnd, badTimes)
                        continue;
                    end

                    allStart(end+1) = s1; %#ok<AGROW>
                    allEnd(end+1)   = s2; %#ok<AGROW>
                    allNote{end+1}  = sprintf('ch=%s', chanNames{iChanList(ii)}); %#ok<AGROW>
                end
            end

            % Ensure events array exists on ORIGINAL raw
            if ~isfield(sFile,'events') || isempty(sFile.events)
                sFile.events = repmat(db_template('event'), 0);
            end

            % Find existing event on ORIGINAL raw
            iEvt = find(strcmpi({sFile.events.label}, evtName), 1);

            % If no detections survived, clear existing output if overwrite
            if isempty(allStart)
                if ~isempty(iEvt) && isOverwrite
                    sEvent = sFile.events(iEvt);
                    sEvent.epochs = [];
                    sEvent.times = [];
                    sEvent.reactTimes = [];
                    sEvent.channels = [];
                    sEvent.notes = [];
                    sFile.events(iEvt) = sEvent;

                    DataMat.F = sFile;
                    bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);
                end

                bst_report('Info', sProcess, sInputs(iFile), 'No HFO candidates detected on filtered signal outside BAD segments.');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            % Build event times [2 x N] USING ORIGINAL RAW TIME AXIS
            evtTimes = [Time(allStart); Time(allEnd)];

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

            % Assign event fields (extended events) ON ORIGINAL RAW
            sEvent.times    = evtTimes;
            sEvent.epochs   = ones(1, size(evtTimes,2));
            sEvent.channels = [];
            sEvent.notes    = allNote;

            sFile.events(iEvt) = sEvent;

            % Save back to ORIGINAL raw file
            DataMat.F = sFile;
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);

            bst_report('Info', sProcess, sInputs(iFile), ...
                sprintf('Detected %d HFO candidates on filtered signal and saved them on the ORIGINAL raw.', size(evtTimes,2)));

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


function badTimes = get_bad_intervals(events)
    badTimes = zeros(2,0);

    if isempty(events)
        return;
    end

    labels = {events.label};
    iBad = find(~cellfun('isempty', regexpi(labels, 'BAD')));

    for i = iBad
        if ~isfield(events(i), 'times') || isempty(events(i).times)
            continue;
        end

        t = events(i).times;

        % Convert simple events to zero-length extended intervals if needed
        if size(t,1) == 1
            t = [t; t];
        elseif size(t,1) ~= 2
            continue;
        end

        t = sort(t, 1);
        badTimes = [badTimes, t]; %#ok<AGROW>
    end

    if ~isempty(badTimes)
        [~, ord] = sort(badTimes(1,:));
        badTimes = badTimes(:, ord);
    end
end


function tf = overlaps_bad(t1, t2, badTimes)
    if isempty(badTimes)
        tf = false;
        return;
    end

    inter = min(t2, badTimes(2,:)) - max(t1, badTimes(1,:));
    tf = any(inter > 0);  % strict overlap only
end