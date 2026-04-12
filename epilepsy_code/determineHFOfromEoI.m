function[HFO] = determineHFOfromEoI(channel,EoI,Fs,time,HFO,time_frequency_map)
    for i = 1:length(EoI)
        if EoI(i)~=0
            %fprintf('\nPosition: %d',i);
            %Compute the average of all the morlet wavelets above 60Hz
            frequencies = 1:1:size(time_frequency_map,1);
            freq_indices = find(frequencies>40);
            high_freq_values = time_frequency_map(freq_indices,:);

            average_morlet = mean(high_freq_values,'all');

            sample_index = i;
            % Define the window around the time point
            window_size = 0.1; %Window size in seconds (e.g 100ms)
            half_window_samples = round(window_size * Fs / 2);
            start_idx = max(1, sample_index - half_window_samples);
            end_idx = min(length(channel), sample_index + half_window_samples);

            % Extract the signal segment
            signal_segment = (channel(start_idx:end_idx));

            % Compute the Power Spectrum using FFT
            N = length(signal_segment); % Number of samples in the segment
            frequencies = (0:N-1) * (Fs / N); % Frequency vector
            spectrum = abs(fft(signal_segment)).^2 / N; % Power spectrum (magnitude squared)
    
            %Convert to decibel units (10 * log10 of power in uV^2/Hz)
            spectrum_db = 10 * log10(spectrum);
            freq_range = frequencies(1:floor(N/2)); % only positive frequencies up to Nyquist
            power_range = spectrum_db(1:floor(N/2)); %Corresponfing power values
            
            %I am gonna put this loop to check if all the power range
            %values are negative (importan for applying the threshold)
            c=0;
            j = 1;
            while j <= length(power_range) && c == 0
                if power_range(j) >= 0
                    c = 1;
                else 
                    c = 0;
                end
                j = j+1;
            end
            if c == 1
                disp('WARNING: not all power values are negative');
            end
            
            %Identify the HiFP (High Frequency peack)
            valid_idx_HiFP = find(40<=freq_range & freq_range<=size(time_frequency_map,1)); %Indices of frequencies >=40 Hz
            [HiFP_value, HiFP_idx_local] = max(power_range(valid_idx_HiFP));% Local index in valid range
            HiFP_idx_global = valid_idx_HiFP(HiFP_idx_local); % Global index in the full range
            HiFP_frequency = freq_range(HiFP_idx_global); % Frequency of the HiFP
            
            %look for the trough
            valid_idx_trough = find(freq_range >= 30 & freq_range <= HiFP_frequency); % Indices between 40 Hz and HiFP
            [trough_value, trough_idx_local] = min(power_range(valid_idx_trough)); % Minimum in this range
            trough_idx_global = valid_idx_trough(trough_idx_local); % Global index
            trough_frequency = freq_range(trough_idx_global); % Frequency of the trough

            %Look for the LoFP (low frequency Peak)
            valid_idx_LoFP = find(freq_range < trough_frequency); % Indices below the trough frequency
            [LoFP_value, LoFP_idx_local] = max(power_range(valid_idx_LoFP)); % Maximum in this range
            LoFP_idx_global = valid_idx_LoFP(LoFP_idx_local); % Global index
            LoFP_frequency = freq_range(LoFP_idx_global); % Frequency of the LoFP
            
            %Determine if the sample is a HFO
            if HiFP_value>0 && trough_value>0
                ratio1 = trough_value/HiFP_value;
            elseif HiFP_value<0 && trough_value<0
                ratio1 = HiFP_value/trough_value;
            elseif HiFP_value>0 && trough_value<0
                ratio1 = abs(trough_value) - abs(HiFP_value);
            end

            if HiFP_value>0 && LoFP_value>0
                ratio2 = HiFP_value/LoFP_value;
            elseif HiFP_value<0 && trough_value<0
                ratio2 = LoFP_value/HiFP_value;
            elseif HiFP_value>0 && trough_value<0
                ratio2 = abs(HiFP_value) - abs(LoFP_value);
            end
                        
            if ratio1<0.95 && ratio2>0.6
                HFO(sample_index) = EoI(sample_index);
            else
                HFO(sample_index) = 0;
            end
            
            
            
            if HFO(sample_index) == 0 && size(time_frequency_map,1)>80
                %Identify the HiFP (High Frequency peack)
                valid_idx_HiFP = find(freq_range>=80 & freq_range<=size(time_frequency_map,1)); %Indices of frequencies >=80 Hz
                [HiFP_value2, HiFP_idx_local2] = max(power_range(valid_idx_HiFP));% Local index in valid range
                HiFP_idx_global2 = valid_idx_HiFP(HiFP_idx_local2); % Global index in the full range
                HiFP_frequency2 = freq_range(HiFP_idx_global2); % Frequency of the HiFP

                %look for the trough
                valid_idx_trough = find(freq_range >= 40 & freq_range <= HiFP_frequency2); % Indices between 40 Hz and HiFP
                [trough_value2, trough_idx_local2] = min(power_range(valid_idx_trough)); % Minimum in this range
                trough_idx_global2 = valid_idx_trough(trough_idx_local2); % Global index
                trough_frequency2 = freq_range(trough_idx_global2); % Frequency of the trough
                
                ratioHFO1 = trough_value2/HiFP_value2;
                ratioHFO2 = LoFP_value/HiFP_value2;
                
                if HiFP_value2 == HiFP_value
                    a = 0;
                else
                    a = 1;
                end
                
                if ratioHFO1>1.1 && ratioHFO2>0.6 && HiFP_value2~=HiFP_value && a==1
                    %disp('here');
                    %fprintf('sample: %d',sample_index);
                    HFO(sample_index) = EoI(sample_index);
                    HiFP_frequency = HiFP_frequency2;
                    if HiFP_value2>HiFP_value
                        HiFP_frequency = HiFP_frequency2;
                        %disp('Entering here');
                        %fprintf('sample: %d',sample_index);
                    end
                else
                    HFO(sample_index) = 0;
                end
            end
    
            n = round(HiFP_frequency);
            morlet_value = time_frequency_map(round(HiFP_frequency),sample_index);

            ratioMorlet = morlet_value/average_morlet;
            max_morlet = max(time_frequency_map(:));
            morlet_threshold = 0.3*max_morlet;
            if HFO(sample_index)~= 0 && ratioMorlet<0.5  
                HFO(sample_index) = 0;
            end

            if HFO(sample_index)~= 0 && morlet_value<morlet_threshold
                HFO(sample_index) = 0;
            end  
        end
    end

%Merge together HFO in an interval less than 30ms (10 samples in between)

%Before merginh HFO let's see if it actualy exists HFO
    exit = 0;
    i = 1;
    while i <= length(HFO) && exit == 0
        if HFO(i) ~= 0
            exit = 1;
        else
            i = i+1;
        end
    end
    
    %Merge HFO
    interval_threshold = (abs(time(1)) - abs(time(12)));
    i = 1;
    while i <= length(HFO) && exit == 1
        if HFO(i)~=0 && i<length(HFO) 
            a = 0;
            %fprintf('\ni before while loop:%d\n',i);
            z = i+1;
            cont = i;
            while a == 0 && z<=length(HFO)
                if HFO(z)~=0
                    z = z + 1;
                    a = 0;
                else
                    a=1;
                    i = z-1;
                    cont = z-1;
                end
            end
            %fprintf('\ni after while loop:%d\n',i);
            %fprintf('cont:%d\n',cont);
            b = cont+1;
            c = 0;
                while b<=length(HFO) && c==0
                    if HFO(b)~=0 
                        if time(b)<0 && time(cont)<0
                            interval = abs(time(cont)) - abs(time(b));
                        else
                            interval = 8000; %rando interval above interval_threshold
                        end
                        %fprintf('For pos %d and pos %d, the interval is: %d\n',cont,b,interval);
                        c = 1;
                        if interval<=interval_threshold
                            for j = (cont+1):b
                                if HFO(j) == 0
                                    HFO(j) = channel(j);
                                    %fprintf('we have include the EoI in pos %d because the interval is %d\n',j,interval);
                                end
                            end
                        end
                    end
                    b = b + 1;
                end
        end
        i = i+1;
    end
    % Discard the HFO that are less than 30ms
    % I decided 30 ms because in the case we are measuring 120 Hz
    % oscillation (120 oscillations in a second), then, 4 oscillations will
    % require 4/120 = 0.03 ms
    i = 1;
    while i <= length(HFO) 
        if HFO(i)~=0
            t1 = time(i);
            a = 0;
            b = i+1; 
            while a == 0 && b<=length(HFO)
                if HFO(b)~=0
                    b=b+1;
                else
                    t2 = time(b);
                    a = 1;
                end
            end
            if a == 0
                t2 = time(b-1);
            end
            HFOduration = t2 - t1;
            if HFOduration < 0.03 && b<=length(HFO)
                for j = i:b
                    HFO(j) = 0;
                end
            end
            i = b+1;
        else
            i = i+1;
        end
    end
    

end
