function[HFO_ascendent_slope,yes_HFO,start_HFO,end_HFO] = determineHFOinSpikeRegion(HFO,time,minimum_index, index_spike)
    HFO_ascendent_slope = zeros(1,length(time));
    for i=minimum_index:index_spike
        HFO_ascendent_slope(i) = HFO(i);
    end
    a = 0;
    if HFO_ascendent_slope(:) == 0
        a = 1;
    end
     %lets see how much time the oscilation last inside the spike
     i = 1;
     k = 0;     
    while i <= length(HFO_ascendent_slope) && a==0
        if HFO_ascendent_slope(i) ~= 0 
            a = 1;
            j = i + 1;
            start_HFO = time(i);
            k = 1;
            m = 0;
            while j <= length(HFO_ascendent_slope) && m == 0
                if HFO_ascendent_slope(j) == 0 
                    end_HFO = time(j-1);
                    m = 1;
                else
                    j = j+1;
                end
            end
        else
            i = i+1;
        end
    end
    if k == 0
        start_HFO = 0;
        end_HFO = 0;
    end
    

        
    yes_HFO = 0;
    i = 1;
    while i<= length(HFO_ascendent_slope) && yes_HFO ==0
        if HFO_ascendent_slope(i)~=0 
            yes_HFO = 1;
        else
            i = i+1;
        end
    end  
    
      
end