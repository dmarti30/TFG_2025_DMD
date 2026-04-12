function varargout =process_evt_detect_overlap_extended( varargin )
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
    sProcess.Comment     = 'Events: Overlap (extended)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 48;

    sProcess.InputTypes  = {'raw','data'};
    sProcess.OutputTypes = {'raw','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Event group A label
    sProcess.options.evta.Comment = 'Event group A (extended label):';
    sProcess.options.evta.Type    = 'text';
    sProcess.options.evta.Value   = '';

    % Event group B label
    sProcess.options.evtb.Comment = 'Event group B (extended label):';
    sProcess.options.evtb.Type    = 'text';
    sProcess.options.evtb.Value   = '';

    % Output label
    sProcess.options.outname.Comment = 'New event label (output):';
    sProcess.options.outname.Type    = 'text';
    sProcess.options.outname.Value   = 'Overlap_AB';

    % Minimum overlap ratio
    sProcess.options.ratio.Comment = 'Min overlap ratio (of shorter event):';
    sProcess.options.ratio.Type    = 'value';
    sProcess.options.ratio.Value   = {0.33, '', 2};

    % Overwrite output group
    sProcess.options.overwrite.Comment = 'Overwrite output label if it already exists';
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

    evtA = strtrim(sProcess.options.evta.Value);
    evtB = strtrim(sProcess.options.evtb.Value);
    outName = strtrim(sProcess.options.outname.Value);
    ratio = sProcess.options.ratio.Value{1};
    isOverwrite = sProcess.options.overwrite.Value;

    if isempty(evtA) || isempty(evtB)
        bst_report('Error', sProcess, [], 'Please enter both event group labels A and B.');
        return;
    end
    if isempty(outName)
        outName = 'Overlap_AB';
    end
    ratio = max(0, min(1, ratio));

    for iFile = 1:numel(sInputs)
        try
            isRaw = strcmpi(sInputs(iFile).FileType, 'raw');

            % ===== Load file and get Events structure =====
            if isRaw
                DataMat = in_bst_data(sInputs(iFile).FileName, 'F');
                sFile = DataMat.F;
                if ~isfield(sFile,'events') || isempty(sFile.events)
                    bst_report('Warning', sProcess, sInputs(iFile), 'No events found in this RAW file.');
                    OutputFiles{iFile} = sInputs(iFile).FileName;
                    continue;
                end
                Events = sFile.events;
            else
                DataMat = in_bst_data(sInputs(iFile).FileName);
                if ~isfield(DataMat,'Events') || isempty(DataMat.Events)
                    bst_report('Warning', sProcess, sInputs(iFile), 'No events found in this DATA file.');
                    OutputFiles{iFile} = sInputs(iFile).FileName;
                    continue;
                end
                Events = DataMat.Events;
            end

            % ===== Find event groups A and B =====
            iA = find(strcmpi({Events.label}, evtA), 1);
            iB = find(strcmpi({Events.label}, evtB), 1);

            if isempty(iA)
                bst_report('Error', sProcess, sInputs(iFile), ['Event group A not found: ' evtA]);
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end
            if isempty(iB)
                bst_report('Error', sProcess, sInputs(iFile), ['Event group B not found: ' evtB]);
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            tA = Events(iA).times;
            tB = Events(iB).times;

            % ===== Enforce EXTENDED events =====
            if isempty(tA) || size(tA,1) ~= 2
                bst_report('Error', sProcess, sInputs(iFile), 'Event group A is not EXTENDED (times must be 2xN).');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end
            if isempty(tB) || size(tB,1) ~= 2
                bst_report('Error', sProcess, sInputs(iFile), 'Event group B is not EXTENDED (times must be 2xN).');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            % Make sure starts <= ends
            tA = sort(tA,1);
            tB = sort(tB,1);

            % Sort by start time
            [~, idxA] = sort(tA(1,:));
            [~, idxB] = sort(tB(1,:));
            tA = tA(:, idxA);
            tB = tB(:, idxB);

            % ===== Compute overlap intervals with ratio criterion =====
            ov = compute_overlap_intervals(tA, tB, ratio);   % 2xK
            if isempty(ov)
                bst_report('Info', sProcess, sInputs(iFile), 'No overlaps meeting the criterion were found.');
                OutputFiles{iFile} = sInputs(iFile).FileName;
                continue;
            end

            % Merge overlaps that strictly overlap (optional de-dup / cleanup)
            ov = merge_intervals_strict(ov);

            % ===== Create / overwrite output event group =====
            iOut = find(strcmpi({Events.label}, outName), 1);

            if isempty(iOut)
                iOut = numel(Events) + 1;
                % Create new event template (works in both RAW and DATA)
                if exist('db_template','file')
                    Events(iOut) = db_template('event');
                else
                    % Minimal fallback structure
                    Events(iOut) = struct('label',[],'color',[],'epochs',[],'times',[], ...
                                          'samples',[],'reactTimes',[],'select',1,'channels',[],'notes',[]);
                end
                Events(iOut).label = outName;

                % Assign color if possible
                if isRaw && exist('panel_record','file')
                    Events(iOut).color = panel_record('GetNewEventColor', iOut, Events);
                else
                    Events(iOut).color = [1 0 0];
                end
            else
                if isOverwrite
                    Events(iOut).times = [];
                    Events(iOut).epochs = [];
                    Events(iOut).reactTimes = [];
                    Events(iOut).samples = [];
                    Events(iOut).channels = [];
                    Events(iOut).notes = [];
                end
            end

            % Append intervals
            Events(iOut).times  = [Events(iOut).times, ov];
            Events(iOut).epochs = [Events(iOut).epochs, ones(1, size(ov,2))];

            % ===== Save back =====
            if isRaw
                sFile.events = Events;
                DataMat.F = sFile;
                bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);
            else
                DataMat.Events = Events;
                bst_save(file_fullpath(sInputs(iFile).FileName), DataMat, 'v6', 1);
            end

            bst_report('Info', sProcess, sInputs(iFile), sprintf('Created %d overlap events in "%s".', size(ov,2), outName));
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


%% ===== LOCAL: compute overlaps with 33% criterion =====
function ov = compute_overlap_intervals(tA, tB, ratio)
% tA, tB: 2xN (sorted by start)
% ratio criterion: overlap >= ratio * min(durA, durB)
% output ov: 2xK with [max(start); min(end)] for each accepted overlap

    i = 1; j = 1;
    nA = size(tA,2);
    nB = size(tB,2);
    ov = zeros(2,0);

    while (i <= nA) && (j <= nB)
        a1 = tA(1,i); a2 = tA(2,i);
        b1 = tB(1,j); b2 = tB(2,j);

        % Compute intersection
        s = max(a1,b1);
        e = min(a2,b2);
        inter = e - s;

        if inter > 0  % IMPORTANT: touching (inter==0) is NOT overlap
            durA = max(eps, a2 - a1);
            durB = max(eps, b2 - b1);
            minDur = min(durA, durB);

            if inter >= (ratio * minDur)
                ov(:,end+1) = [s; e]; %#ok<AGROW>
            end
        end

        % Advance pointer
        if a2 < b2
            i = i + 1;
        else
            j = j + 1;
        end
    end
end


%% ===== LOCAL: merge overlap intervals (strict overlap only) =====
function out = merge_intervals_strict(ov)
% ov: 2xK, unsorted ok. Merge only if next_start < current_end (strict overlap).
    if isempty(ov)
        out = ov;
        return;
    end
    [~, idx] = sort(ov(1,:));
    ov = ov(:,idx);

    curS = ov(1,1);
    curE = ov(2,1);

    out = zeros(2,0);

    for k = 2:size(ov,2)
        s = ov(1,k);
        e = ov(2,k);

        if s < curE  % strict overlap (do NOT merge if s == curE)
            curE = max(curE, e);
        else
            out(:,end+1) = [curS; curE]; %#ok<AGROW>
            curS = s; curE = e;
        end
    end

    out(:,end+1) = [curS; curE];
end
