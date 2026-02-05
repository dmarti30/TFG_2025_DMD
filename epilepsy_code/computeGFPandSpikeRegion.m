function[GFP,minimum_index, index_spike, min_after_max] = computeGFPandSpikeRegion(FNormal,timeNormal)
%     Fnew = zeros(size(FNormal,1)-2,size(FNormal,2));
%         for i = 1:(size(FNormal,1)-2)
%             for j = 1:size(FNormal,2)
%                 Fnew(i,j) = FNormal(i,j);
%             end
%         end

        GFP = std(FNormal,1);
            %Let's look for the samples corresponding to the ascendent slope of the
        %spike

        % Look for the sample where t = 0
        [~, index_t0] = min(abs(timeNormal - 0)); % sample_t0 = closest sample to t = 0

        % Find all the local minimums of the GFP
        local_min_index = islocalmin(GFP(1:index_t0)); %return a vector with 0 and 1 when a min is found
        
        %Look for the min just before the zero
        for i = 1:length(local_min_index)
            if local_min_index(i) == 1
                minimum_index = i;
            end
        end

  
        %Look for the maximum after minimum_index
        index_spike_array = islocalmax(GFP(minimum_index:index_t0+10)); %index_spike will be the local index between this range
        i = 2;
        value = GFP(index_t0);
        a = 0;
        while i<=length(index_spike_array) 
            if index_spike_array(i) == 1 && GFP(minimum_index + (i-1))>value %to make sure we found the bigger maximum between all the local maximumns within that range
                index_spike = i;
                value = GFP(minimum_index + (index_spike -1));
                a = 1;
            else
                i = i+1;
            end
        end
        if a == 0
            index_spike = index_t0;
        end
        if a == 1
            index_spike = minimum_index + (index_spike -1); %to find the corresponding index along the whole range
        end

        if abs(timeNormal(index_spike) - timeNormal(minimum_index))<0.025
            i = minimum_index;
            a = 0;
            while i>0 && a==0
                if local_min_index(i) == 1 && abs(timeNormal(index_spike) - timeNormal(i))>0.025 && GFP(i)<GFP(minimum_index) && abs(timeNormal(index_spike) - timeNormal(i))<0.07
                    minimum_index = i;
                    a = 1;
                else
                    i = i-1;
                end
            end
        end
        
        min_after_max_array = islocalmin(GFP(index_spike:length(timeNormal)));
        i = 1;
        a = 0;
        while i <= length(min_after_max_array) && a==0
            if min_after_max_array(i) == 1
                min_after_max = i + index_spike - 1;
                if abs(timeNormal(min_after_max) - timeNormal(index_spike))>0.025 && GFP(i)<GFP(min_after_max) && abs(timeNormal(min_after_max) - timeNormal(index_spike))<0.07
                    a = 1;
                else
                    i = i+1;
                end
            else
                i = i+1;
            end
        end   
end