function varargout = process_evt_detect_overlap_keep_spike(varargin)
% PROCESS_EVT_DETECT_OVERLAP_KEEP_SPIKE:
% Detect temporal overlap between two extended event groups (typically IEDs and HFOs),
% but instead of saving only the overlapping portion, create a new event group that
% keeps the FULL spike/IED interval whenever that spike contains an HFO with sufficient overlap.
%
% Overlap criterion:
%   overlap_duration >= min_overlap_pct * min(spike_duration, hfo_duration)
%
% Example:
%   Spike interval: [0.900 1.000]
%   HFO interval:   [0.960 0.985]
%   If overlap passes the threshold, output event = [0.900 1.000]
%
% This function is designed for Brainstorm custom processes.

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Overlap analysis (keep full spike)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 59;

    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Spike / IED event label
    sProcess.options.spikelabel.Comment = 'Spike / IED event label:';
    sProcess.options.spikelabel.Type    = 'text';
    sProcess.options.spikelabel.Value   = 'IED';

    % HFO event label
    sProcess.options.hfolabel.Comment = 'HFO event label:';
    sProcess.options.hfolabel.Type    = 'text';
    sProcess.options.hfolabel.Value   = 'HFO';

    % New output label
    sProcess.options.newlabel.Comment = 'Output event label:';
    sProcess.options.newlabel.Type    = 'text';
    sProcess.options.newlabel.Value   = 'Spike+HFO';

    % Minimum overlap percentage
    sProcess.options.minoverlap.Comment = 'Minimum overlap (% of shorter event):';
    sProcess.options.minoverlap.Type    = 'value';
    sProcess.options.minoverlap.Value   = {33, '%', 0};

    % Overwrite existing output event
    sProcess.options.overwrite.Comment = 'Overwrite existing output event';
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

    spikeLabel   = strtrim(sProcess.options.spikelabel.Value);
    hfoLabel     = strtrim(sProcess.options.hfolabel.Value);
    outLabel     = strtrim(sProcess.options.newlabel.Value);
    minOverlapPc = sProcess.options.minoverlap.Value{1};
    isOverwrite  = sProcess.options.overwrite.Value;

    if isempty(spikeLabel)
        spikeLabel = 'IED';
    end
    if isempty(hfoLabel)
        hfoLabel = 'HFO';
    end
    if isempty(outLabel)
        outLabel = 'Spike+HFO';
    end
    if isempty(minOverlapPc) || ~isnumeric(minOverlapPc)
        minOverlapPc = 33;
    end
    minOverlapFrac = max(0, minOverlapPc) / 100;

    for iFile = 1:numel(sInputs)
        try
            DataMat = in_bst_data(sInputs(iFile).FileName, 'F');
            sFile   = DataMat.F;

            if ~isfield(sFile, 'events') || isempty(sFile.events)
                bst_report('Error', sProcess, sInputs(iFile), 'This RAW file does not contain any events.');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            iSpike = find_event_by_label(sFile.events, spikeLabel);
            iHFO   = find_event_by_label(sFile.events, hfoLabel);

            if isempty(iSpike)
                bst_report('Error', sProcess, sInputs(iFile), sprintf('Event group "%s" not found.', spikeLabel));
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end
            if isempty(iHFO)
                bst_report('Error', sProcess, sInputs(iFile), sprintf('Event group "%s" not found.', hfoLabel));
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            spikeTimes = normalize_event_times(sFile.events(iSpike).times);
            hfoTimes   = normalize_event_times(sFile.events(iHFO).times);

            nSpike = size(spikeTimes, 2);
            nHFO   = size(hfoTimes, 2);

            if nSpike == 0
                bst_report('Warning', sProcess, sInputs(iFile), sprintf('Event group "%s" is empty.', spikeLabel));
            end
            if nHFO == 0
                bst_report('Warning', sProcess, sInputs(iFile), sprintf('Event group "%s" is empty.', hfoLabel));
            end

            keepSpike = false(1, nSpike);

            for iSpk = 1:nSpike
                s1 = spikeTimes(1, iSpk);
                s2 = spikeTimes(2, iSpk);
                spikeDur = max(eps, s2 - s1);

                for iEvt = 1:nHFO
                    h1 = hfoTimes(1, iEvt);
                    h2 = hfoTimes(2, iEvt);
                    hfoDur = max(eps, h2 - h1);

                    overlapDur = min(s2, h2) - max(s1, h1);

                    if overlapDur > 0
                        minRequired = minOverlapFrac * min(spikeDur, hfoDur);
                        if overlapDur >= minRequired
                            keepSpike(iSpk) = true;
                            break;
                        end
                    end
                end
            end

            outTimes = spikeTimes(:, keepSpike);
            outTimes = unique_event_columns(outTimes);

            % Ensure events exist
            if ~isfield(sFile, 'events') || isempty(sFile.events)
                sFile.events = repmat(db_template('event'), 0);
            end

            % Find/create output event
            iOut = find_event_by_label(sFile.events, outLabel);

            if isempty(iOut)
                iOut = numel(sFile.events) + 1;
                sEvent = db_template('event');
                sEvent.label = outLabel;
                sEvent.color = panel_record('GetNewEventColor', iOut, sFile.events);
            else
                sEvent = sFile.events(iOut);
            end

            if isOverwrite
                sEvent.times      = [];
                sEvent.epochs     = [];
                sEvent.reactTimes = [];
                sEvent.channels   = [];
                sEvent.notes      = [];
            end

            % Append full spike intervals that matched at least one HFO
            if ~isempty(outTimes)
                sEvent.times  = [normalize_event_times(sEvent.times), outTimes];
                sEvent.times  = unique_event_columns(sEvent.times);
                sEvent.epochs = ones(1, size(sEvent.times, 2));
            else
                % Keep empty output group visible in the event list
                sEvent.times  = normalize_event_times(sEvent.times);
                sEvent.epochs = ones(1, size(sEvent.times, 2));
            end

            sEvent.channels   = [];
            sEvent.reactTimes = [];
            sEvent.notes      = [];

            sFile.events(iOut) = sEvent;

            % Save back to RAW file
            DataMat.F = sFile;
            bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);

            bst_report('Info', sProcess, sInputs(iFile), ...
                sprintf('Created "%s" with %d spike(s) containing HFO overlap.', outLabel, size(outTimes, 2)));

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
function iEvt = find_event_by_label(events, label)
    iEvt = find(strcmpi({events.label}, label), 1);
end


function times = normalize_event_times(timesIn)
    if isempty(timesIn)
        times = zeros(2,0);
        return;
    end

    times = timesIn;

    if size(times,1) == 1
        % Single events -> convert to zero-duration extended events
        times = [times; times];
    elseif size(times,1) == 2
        % OK
    elseif size(times,2) == 2
        % Nx2 -> convert to 2xN
        times = times';
    else
        error('Unsupported event time format.');
    end

    times = sort(times, 1);
end


function times = unique_event_columns(timesIn)
    times = normalize_event_times(timesIn);

    if isempty(times)
        return;
    end

    % Remove exact duplicate columns while preserving order
    keep = true(1, size(times,2));
    seen = containers.Map();

    for i = 1:size(times,2)
        key = sprintf('%.12f_%.12f', times(1,i), times(2,i));
        if isKey(seen, key)
            keep(i) = false;
        else
            seen(key) = true;
        end
    end

    times = times(:, keep);
end