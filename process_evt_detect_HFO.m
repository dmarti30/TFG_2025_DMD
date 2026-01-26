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
    sProcess.Comment     = 'Detect HFO in Spikes';
    sProcess.FileTag     = 'HFO';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'HFO detection';
    sProcess.Index       = 55;
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % Configurable parameters, example:
    sProcess.options.Fs.Comment = 'Sampling frequency (Hz)';
    sProcess.options.Fs.Type    = 'value';
    sProcess.options.Fs.Value   = 512; % Default value
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    % Extract parameters
    Fs = sProcess.options.Fs.Value;
    nFiles = length(sInputs);
    
    % Initialize variables for global results
    all_HFO = [];
    all_HFO_spike_region = [];
    all_EoI = [];
    all_GFP = [];
    all_start_end_spike_region = [];
    all_start_end_HFO_within_spike = [];
    spikes_with_HFO_spike_region = [];
    channels_HFO = [];
    
    % Process each file (spike)
    for f = 1:nFiles
        % Load data
        inMat = in_bst(sInputs(f).FileName);
        
        % Assuming typical structure, modify if necessary
        F = inMat.F; % Channels x Samples
        time = inMat.Time;
        
        if isfield(inMat, 'Morlets')
            morlets = inMat.Morlets; % Channels x Samples x Frequencies
        else
            error('Morlets matrix missing in input file.');
        end
        
        nChannels = size(F,1);
        nSamples = size(F,2);
        nFreq = size(morlets,3);
        
        % Initialize for this spike
        EoI = zeros(nChannels,nSamples);
        HFO = zeros(nChannels,nSamples);
        HFO_spike_region = zeros(nChannels,nSamples);
        
        % Detect Events of Interest (EoI)
        for ch = 1:nChannels
            EoI(ch,:) = determineEoI(F(ch,:), time);
        end
        
        % Determine if EoI is HFO with Morlets
        for ch = 1:nChannels
            tf_map = squeeze(morlets(ch,:,:))';
            HFOrow = zeros(1,nSamples);
            HFOrow = determineHFOfromEoI(F(ch,:), EoI(ch,:), Fs, time, HFOrow, tf_map);
            HFO(ch,:) = HFOrow;
        end
        
        % Calculate GFP and spike region
        [GFP, minimum_index, spike_index, min_after_spike] = computeGFPandSpikeRegion(F, time);
        
        % Determine HFO in spike region
        channels_HFO_spike = zeros(1,nChannels);
        all_start_end_HFO_within_spike_current = zeros(nChannels,2);
        for ch = 1:nChannels
            [HFO_spike_row, yes_HFO, start_HFO, end_HFO] = determineHFOinSpikeRegion(HFO(ch,:), time, minimum_index, min_after_spike);
            HFO_spike_region(ch,:) = HFO_spike_row;
            channels_HFO_spike(ch) = yes_HFO;
            all_start_end_HFO_within_spike_current(ch,1) = start_HFO;
            all_start_end_HFO_within_spike_current(ch,2) = end_HFO;
        end
        
        % Save global results
        if isempty(all_HFO)
            all_HFO = zeros(nFiles, nChannels, nSamples);
            all_EoI = zeros(nFiles, nChannels, nSamples);
            all_HFO_spike_region = zeros(nFiles, nChannels, nSamples);
            all_GFP = zeros(nFiles, nSamples);
            all_start_end_spike_region = zeros(nFiles,3);
            all_start_end_HFO_within_spike = zeros(nFiles, nChannels, 2);
            spikes_with_HFO_spike_region = zeros(nFiles,1);
            channels_HFO = zeros(nFiles, nChannels);
        end
        
        all_HFO(f,:,:) = HFO;
        all_EoI(f,:,:) = EoI;
        all_HFO_spike_region(f,:,:) = HFO_spike_region;
        all_GFP(f,:) = GFP;
        all_start_end_spike_region(f,:) = [minimum_index, min_after_spike, spike_index];
        all_start_end_HFO_within_spike(f,:,:) = all_start_end_HFO_within_spike_current;
        channels_HFO(f,:) = channels_HFO_spike;
        spikes_with_HFO_spike_region(f) = any(channels_HFO_spike);
    end
    
    % Post-processing: remove invalid HFO (start==end or too close)
    for i = 1:length(spikes_with_HFO_spike_region)
        if spikes_with_HFO_spike_region(i)
            for ch = 1:nChannels
                if channels_HFO(i,ch) == 1
                    start_hfo = all_start_end_HFO_within_spike(i,ch,1);
                    end_hfo = all_start_end_HFO_within_spike(i,ch,2);
                    if abs(start_hfo - end_hfo) < 0.01
                        channels_HFO(i,ch) = 0;
                    end
                end
            end
            if all(channels_HFO(i,:) == 0)
                spikes_with_HFO_spike_region(i) = 0;
            end
        end
    end
    
    % Save file with results
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sInputs(1).FileName), 'HFO_detected.mat');
    save(OutputFile, 'all_HFO', 'all_HFO_spike_region', 'all_EoI', 'all_GFP', ...
        'all_start_end_spike_region', 'all_start_end_HFO_within_spike', 'spikes_with_HFO_spike_region', 'channels_HFO', '-v7.3');
    
    OutputFiles = {OutputFile};
end
